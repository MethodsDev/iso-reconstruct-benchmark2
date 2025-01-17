version 1.0

# This task uses ESPRESSO version 1.5.0
task espressoTask {
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
        Int diskSizeGB = 256
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/espresso:latest"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String OutDir = "ESPRESSO_out"
    String samples_filename = "espresso_samples.tsv"

    command <<<
        bash ~{monitoringScript} > monitoring.log &
        mkdir -p ESPRESSO_out
        mkdir -p ESPRESSO_out/ID
        mkdir -p ESPRESSO_out/Quant

        samtools view -h -o input.sam ~{inputBAM}
        
        echo -e "input.sam\tespresso" > ESPRESSO_out/ID/~{samples_filename}
        mkdir -p ESPRESSO_out/ID
        echo -e "input.sam\tespresso" > ESPRESSO_out/ID/~{samples_filename}

        echo -e "input.sam\tespresso" > ESPRESSO_out/Quant/~{samples_filename}
        mkdir -p ESPRESSO_out/Quant
        echo -e "input.sam\tespresso" > ESPRESSO_out/Quant/~{samples_filename}
        
        if [[ "~{referenceAnnotation_reduced}" != "" && ("~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both") ]]; then
            perl /opt/conda/envs/espresso_env/bin/ESPRESSO_S.pl --sort_buffer_size ~{memoryGB} -L ESPRESSO_out/ID/~{samples_filename} -F ~{referenceGenome} -A ~{referenceAnnotation_reduced} -O ESPRESSO_out/ID -T ~{numThreads}
            perl /opt/conda/envs/espresso_env/bin/ESPRESSO_C.pl --sort_buffer_size ~{memoryGB} -I ESPRESSO_out/ID -F ~{referenceGenome} -X 0 -T ~{numThreads}
            perl /opt/conda/envs/espresso_env/bin/ESPRESSO_Q.pl -L ESPRESSO_out/ID/espresso_samples.tsv.updated -A ~{referenceAnnotation_reduced} -T ~{numThreads}            
            mv ESPRESSO_out/ID/espresso_samples_N2_R0_updated.gtf ESPRESSO_out/ID/espressoReducedGTF.gtf
            mv ESPRESSO_out/ID/espresso_samples_N2_R0_abundance.esp ESPRESSO_out/ID/espressoReducedGTFCounts.txt
        fi


        if [[("~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both")]]; then
            perl /opt/conda/envs/espresso_env/bin/ESPRESSO_S.pl --sort_buffer_size ~{memoryGB} -L ESPRESSO_out/Quant/~{samples_filename} -F ~{referenceGenome} -A ~{referenceAnnotation_full} -O ESPRESSO_out/Quant -T ~{numThreads}
            perl /opt/conda/envs/espresso_env/bin/ESPRESSO_C.pl --sort_buffer_size ~{memoryGB} -I ESPRESSO_out/Quant -F ~{referenceGenome} -X 0 -T ~{numThreads}
            perl /opt/conda/envs/espresso_env/bin/ESPRESSO_Q.pl -L ESPRESSO_out/Quant/espresso_samples.tsv.updated -A ~{referenceAnnotation_full} -T ~{numThreads}            
            mv ESPRESSO_out/Quant/espresso_samples_N2_R0_updated.gtf ESPRESSO_out/Quant/espressoFullGTF.gtf
            mv ESPRESSO_out/Quant/espresso_samples_N2_R0_abundance.esp ESPRESSO_out/Quant/espressoCounts.txt
        fi

    >>>

    output {
        File? espressoReducedGTFCounts = "~{OutDir}/ID/espressoReducedGTFCounts.txt"
        File? espressoReducedGTF = "~{OutDir}/ID/espressoReducedGTF.gtf"
        File? espressoCounts = "~{OutDir}/Quant/espressoCounts.txt"
        File? espressoFullGTF = "~{OutDir}/Quant/espressoFullGTF.gtf"

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

workflow espressoWorkflow {
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

    call espressoTask {
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
        File? espressoReducedGTFCounts = espressoTask.espressoReducedGTFCounts
        File? espressoReducedGTF = espressoTask.espressoReducedGTF
        File monitoringLog = espressoTask.monitoringLog
        File? espressoCounts = espressoTask.espressoCounts
        File? espressoFullGTF = espressoTask.espressoFullGTF
    }
}
