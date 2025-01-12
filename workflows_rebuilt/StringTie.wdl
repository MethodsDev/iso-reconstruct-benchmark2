version 1.0

# This task uses StringTie version 2.2.1
task stringtieTask {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String dataType
        String ID_or_Quant_or_Both
        String Reffree_or_Refguided_or_Both
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 250
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/stringtie@sha256:fb579fc315d6976ccfa4094e7b0e5fe45587f56b02a6124bcce463023ead7d5d"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }
    
    String OutDir = "StringTie_out"

    command <<<
        bash ~{monitoringScript} > monitoring.log &
        mkdir -p ~{OutDir}
        mkdir -p ~{OutDir}_IDreffreeQuant
        mkdir -p ~{OutDir}_IDrefguidedQuant

        if [[ "~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            if [[ "~{Reffree_or_Refguided_or_Both}" == "Reffree" || "~{Reffree_or_Refguided_or_Both}" == "Both" ]]; then
                stringtie \
                -o "~{OutDir}/stringtieGTF.gtf" \
                -p ~{numThreads} \
                -L ~{inputBAM} \
                --ref ~{referenceGenome}

                stringtie -e -o "~{OutDir}_IDreffreeQuant/StringTie_quant.gtf" -G "~{OutDir}/stringtieGTF.gtf" -p ~{numThreads} -L ~{inputBAM} --ref ~{referenceGenome}
                echo -e "stringtie\t~{OutDir}_IDreffreeQuant/StringTie_quant.gtf" > ~{OutDir}/stringtie_sample_list2.txt
                prepDE.py -i ~{OutDir}/stringtie_sample_list2.txt -g ~{OutDir}_IDreffreeQuant/gene_count_matrix.csv -t ~{OutDir}/stringtieGTFCounts.csv
            fi

            if [[ "~{Reffree_or_Refguided_or_Both}" == "Refguided" || "~{Reffree_or_Refguided_or_Both}" == "Both" ]]; then
                if [[ -n "~{referenceAnnotation_reduced}" ]]; then
                    stringtie \
                    -o "~{OutDir}/stringtieReducedGTF.gtf" \
                    -G ~{referenceAnnotation_reduced} \
                    -p ~{numThreads} \
                    -L ~{inputBAM} \
                    --ref ~{referenceGenome}

                    stringtie -e -o "~{OutDir}_IDrefguidedQuant/StringTie_quant.gtf" -G "~{OutDir}/stringtieReducedGTF.gtf" -p ~{numThreads} -L ~{inputBAM} --ref ~{referenceGenome}
                    echo -e "stringtie\t~{OutDir}_IDrefguidedQuant/StringTie_quant.gtf" > ~{OutDir}/stringtie_sample_list3.txt
                    prepDE.py -i ~{OutDir}/stringtie_sample_list3.txt -g ~{OutDir}_IDrefguidedQuant/gene_count_matrix.csv -t ~{OutDir}/stringtieReducedGTFCounts.csv
                fi
            fi
        fi

        if [[ "~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            if [[ -n "~{referenceAnnotation_full}" ]]; then
                stringtie -e -o "~{OutDir}/StringTie_quant.gtf" -G ~{referenceAnnotation_full} -p ~{numThreads} -L ~{inputBAM} --ref ~{referenceGenome}
                echo -e "stringtie\t~{OutDir}/StringTie_quant.gtf" > ~{OutDir}/stringtie_sample_list.txt
                prepDE.py -i ~{OutDir}/stringtie_sample_list.txt -g ~{OutDir}/gene_count_matrix.csv -t ~{OutDir}/stringtieCounts.csv
            fi
        fi
    >>>
    
    output {
        File? stringtieGTF = "~{OutDir}/stringtieGTF.gtf"
        File? stringtieReducedGTF = "~{OutDir}/stringtieReducedGTF.gtf"
        File? stringtieCounts = "~{OutDir}/stringtieCounts.csv"
        File monitoringLog = "monitoring.log"
        File? stringtieGTFCounts = "~{OutDir}/stringtieGTFCounts.csv"
        File? stringtieReducedGTFCounts = "~{OutDir}/stringtieReducedGTFCounts.csv"
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
        String Reffree_or_Refguided_or_Both
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
            ID_or_Quant_or_Both = ID_or_Quant_or_Both,
            Reffree_or_Refguided_or_Both = Reffree_or_Refguided_or_Both
    }

    output {
        File? stringtieGTF = stringtieTask.stringtieGTF
        File? stringtieReducedGTF = stringtieTask.stringtieReducedGTF
        File? stringtieCounts = stringtieTask.stringtieCounts
        File monitoringLog = stringtieTask.monitoringLog
        File? stringtieGTFCounts = stringtieTask.stringtieGTFCounts
        File? stringtieReducedGTFCounts = stringtieTask.stringtieReducedGTFCounts
    }
}
