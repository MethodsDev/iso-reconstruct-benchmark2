version 1.0

# This task uses Isosceles version
task isoscelesTask {
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
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isosceles:latest"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }
    
    
    String OutDir = "Isosceles_out"

    command <<<
        bash ~{monitoringScript} > monitoring.log &
        mkdir -p $OutDir

        if [[ "~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
#            stringtie \
#            -o "~{OutDir}/StringTie.gtf" \
#            -p ~{numThreads} \
#            -L ~{inputBAM} \
#            --ref ~{referenceGenome}
            
            if [[ -n "~{referenceAnnotation_reduced}" ]]; then
                stringtie \
                -o "~{OutDir}/StringTie_reduced.gtf" \
                -G ~{referenceAnnotation_reduced} \
                -p ~{numThreads} \
                -L ~{inputBAM} \
                --ref ~{referenceGenome}

            fi
        fi

        if [[ "~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            if [[ -n "~{referenceAnnotation_full}" ]]; then
                stringtie -e -o "~{OutDir}/StringTie_quant.gtf" -G ~{referenceAnnotation_full} -p ~{numThreads} -L ~{inputBAM} --ref ~{referenceGenome}
                echo -e "stringtie\t~{OutDir}/StringTie_quant.gtf" > ~{OutDir}/stringtie_sample_list.txt
                prepDE.py -i ~{OutDir}/stringtie_sample_list.txt -g ~{OutDir}/gene_count_matrix.csv -t ~{OutDir}/StringTie_quant.csv
            fi
        fi
    >>>
    
    output {
        File? stringtieGTF = "~{OutDir}/StringTie.gtf"
        File? stringtieReducedGTF = "~{OutDir}/StringTie_reduced.gtf"
        File? stringtieCounts = "~{OutDir}/StringTie_quant.csv"
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

workflow stringtieWorkflow {
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

    call stringtieTask {
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
        File? stringtieGTF = stringtieTask.stringtieGTF
        File? stringtieReducedGTF = stringtieTask.stringtieReducedGTF
        File? stringtieCounts = stringtieTask.stringtieCounts
        File monitoringLog = stringtieTask.monitoringLog
    }
}
