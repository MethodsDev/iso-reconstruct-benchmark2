version 1.0

# This task uses Oarfish version 4.5.0
task oarfishTask {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String dataType
        String ID_or_Quant_or_Both
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 250
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/oarfish@sha256:84b87e20aa10e2267e117da7b02c8dfa72f23eba6198ca105d26025dece0058a"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }
    String OutDir = "Oarfish_out"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        if [[ "~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            if [[ -z "~{referenceAnnotation_full}" ]]; then
                echo "Full annotation is not provided. Please provide the full annotation."
                exit 1
            fi  
            mkdir -p ~{OutDir}

            samtools bam2fq ~{inputBAM} > ~{OutDir}/tmp.fastq

            gffread ~{referenceAnnotation_full} -g ~{referenceGenome} -w  ~{OutDir}/transcriptome.fa 

            minimap2 -a -t ~{numThreads} ~{OutDir}/transcriptome.fa  ~{OutDir}/tmp.fastq | samtools view --threads ~{numThreads} -bS | samtools sort --threads ~{numThreads} > ~{OutDir}/mapped.sorted.bam
            samtools index ~{OutDir}/mapped.sorted.bam
            
            cd ~{OutDir}
            oarfish --alignments mapped.sorted.bam --output ~{OutDir} --model-coverage --threads ~{numThreads}
            
            mv Oarfish_out.quant Oarfish_quant.tsv
        fi
    >>>

    output {
        File? oarfishCounts = "~{OutDir}/Oarfish_quant.tsv"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}


workflow oarfishWorkflow {
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

    call oarfishTask {
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
        File? oarfishCounts = oarfishTask.oarfishCounts
        File monitoringLog = oarfishTask.monitoringLog
    }
}
