version 1.0

# This task uses TALON version 6.0
task talonTask {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String dataType
        String ID_or_Quant_or_Both
        Int cpu = 8
        Int numThreads = 16
        Int memoryGB = 128
        Int diskSizeGB = 1024
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/talon@sha256:a1d019708bd6b98c456011e56f24d9799aeb6ad45c3c5a33bab931d396b9c0b2"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String OutDir = "TALON_out"
    String talonPrefix = "TALON"
    String talonPrefix2 = "TALON2"

    String datasetName = "TALON_dataset"
    String datasetName2 = "TALON_dataset2"

    command <<<
        bash ~{monitoringScript} > monitoring.log &
        mkdir -p ~{OutDir}
        cd ~{OutDir}
        talon_label_reads --f ~{inputBAM} --t 1 --o ~{talonPrefix} --g ~{referenceGenome}
        samtools calmd -@ ~{numThreads} --reference ~{referenceGenome} "~{talonPrefix}_labeled.sam" > "~{talonPrefix}_labeled.md.sam"

        if [[ "~{ID_or_Quant_or_Both}" == "ID" ]] || [[ "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            if [[ "~{referenceAnnotation_reduced}" != "" ]]; then

                talon_initialize_database --f ~{referenceAnnotation_reduced} --g ~{datasetName} --a ~{datasetName} --o ~{datasetName}
                echo ~{datasetName},~{datasetName},~{dataType},"~{talonPrefix}_labeled.md.sam" > "~{talonPrefix}.csv"
                talon --build ~{datasetName} --db "~{datasetName}.db" --o "~{talonPrefix}_raw" --f "~{talonPrefix}.csv" --threads ~{numThreads}
                talon_filter_transcripts --db "~{datasetName}.db" -a ~{datasetName} --datasets ~{datasetName} --o "~{talonPrefix}_filter"
                talon_create_GTF --build ~{datasetName} --db "~{datasetName}.db" -a ~{datasetName} --o ~{talonPrefix} --whitelist "~{talonPrefix}_filter"
                talon_abundance --db "~{datasetName}.db" --whitelist "~{talonPrefix}_filter" --o ~{talonPrefix}_Quant --build ~{datasetName} -a ~{datasetName}
                ls -l
                mv TALON_talon.gtf talonReducedGTF.gtf
                mv TALON_Quant_talon_abundance_filtered.tsv talonReducedGTFCounts.tsv
            fi
        fi


        if [[ "~{ID_or_Quant_or_Both}" == "Quant" ]] || [[ "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            if [[ "~{referenceAnnotation_full}" != "" ]]; then

                talon_initialize_database --f ~{referenceAnnotation_full} --g ~{datasetName2} --a ~{datasetName2} --o ~{datasetName2}
                echo ~{datasetName2},~{datasetName2},~{dataType},"~{talonPrefix}_labeled.md.sam" > "~{talonPrefix2}.csv"
                talon --build ~{datasetName2} --db "~{datasetName2}.db" --o "~{talonPrefix2}_raw" --f "~{talonPrefix2}.csv" --threads ~{numThreads}
                talon_filter_transcripts --db "~{datasetName2}.db" -a ~{datasetName2} --datasets ~{datasetName2} --o "~{talonPrefix2}_filter"
                talon_create_GTF --build ~{datasetName2} --db "~{datasetName2}.db" -a ~{datasetName2} --o ~{talonPrefix2} --whitelist "~{talonPrefix2}_filter"
                talon_abundance --db "~{datasetName2}.db" --whitelist "~{talonPrefix2}_filter" --o ~{talonPrefix2}_Quant --build ~{datasetName2} -a ~{datasetName2}
                ls -l
                mv TALON2_talon.gtf talonFullGTF.gtf
                mv TALON2_Quant_talon_abundance_filtered.tsv talonCounts.tsv
            fi
        fi

    >>>

    output {
        File? talonReducedGTF = "TALON_out/talonReducedGTF.gtf"
        File monitoringLog = "monitoring.log"
        File? talonReducedGTFCounts = "TALON_out/talonReducedGTFCounts.tsv"
        File? talonCounts = "TALON_out/talonCounts.tsv"
        File? talonFullGTF = "TALON_out/talonFullGTF.gtf"

    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}

workflow talonWorkflow {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String dataType
        String ID_or_Quant_or_Both
    }

    call talonTask {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }

    output {
        File? talonReducedGTF = talonTask.talonReducedGTF
        File monitoringLog = talonTask.monitoringLog
        File? talonReducedGTFCounts = talonTask.talonReducedGTFCounts
        File? talonCounts = talonTask.talonCounts
        File? talonFullGTF = talonTask.talonFullGTF
    }
}
