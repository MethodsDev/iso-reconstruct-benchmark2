version 1.0

# This task uses ESPRESSO version 1.4.0
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
        Int diskSizeGB = 250
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/espresso@sha256:f538303f6457c55e7b3c2a45081e6d8e3053e6f76e56bc65631b7f4aa290b026"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String OutDir = "ESPRESSO_out"
    String samples_filename = "espresso_samples.tsv"

    command <<<
        bash ~{monitoringScript} > monitoring.log &
        mkdir -p ESPRESSO_out
        
        # Convert BAM to SAM
        samtools view -h -o input.sam ~{inputBAM}
        
        # Create espresso_samples.tsv
        echo -e "input.sam\tespresso" > ESPRESSO_out/ID/~{samples_filename}
        
        if [[ "~{referenceAnnotation_reduced}" != "" && ("~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both" || "~{ID_or_Quant_or_Both}" == "ID") ]]; then
            mkdir -p ESPRESSO_out/ID
            perl /usr/src/app/espresso/src/ESPRESSO_S.pl --sort_buffer_size ~{memoryGB} -L ESPRESSO_out/ID/~{samples_filename} -F ~{referenceGenome} -A ~{referenceAnnotation_reduced} -O ESPRESSO_out/ID -T ~{numThreads}
            perl /usr/src/app/espresso/src/ESPRESSO_C.pl --sort_buffer_size ~{memoryGB} -I ESPRESSO_out/ID -F ~{referenceGenome} -X 0 -T ~{numThreads}
            perl /usr/src/app/espresso/src/ESPRESSO_Q.pl -L ESPRESSO_out/ID/espresso_samples.tsv.updated -A ~{referenceAnnotation_reduced} -T ~{numThreads}            
            mv ESPRESSO_out/ID/espresso_samples_N2_R0_updated.gtf ESPRESSO_out/ID/ESPRESSO_reduced.gtf
        fi
    
 #       if [[ "~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
 #           mkdir -p ~{OutDir}/Quant
 #           perl /usr/src/app/espresso/src/ESPRESSO_S.pl --sort_buffer_size ~{memoryGB} -L ~{OutDir}/Quant/~{samples_filename} -F ~{referenceGenome} -A ~{referenceAnnotation_full} -O ~{OutDir}/Quant -T ~{numThreads}
 #           perl /usr/src/app/espresso/src/ESPRESSO_C.pl --sort_buffer_size ~{memoryGB} -I ~{OutDir}/Quant -F ~{referenceGenome} -X 0 -T ~{numThreads}
 #           perl /usr/src/app/espresso/src/ESPRESSO_Q.pl -L ~{OutDir}/Quant/espresso_samples.tsv.updated -A ~{referenceAnnotation_full} -T ~{numThreads}
 #           mv ~{OutDir}/Quant/espresso_samples_N2_R0_abundance.esp ~{OutDir}/Quant/ESPRESSO_quant.txt
 #       fi
    >>>

    output {
        File? espressoCounts = "~{OutDir}/Quant/ESPRESSO_quant.txt"
        File? espressoReducedGTF = "~{OutDir}/ID/ESPRESSO_reduced.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: 16
        memory: "256 GiB"
        disks: "local-disk 500 HDD"
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
        File? espressoCounts = espressoTask.espressoCounts
        File? espressoReducedGTF = espressoTask.espressoReducedGTF
        File monitoringLog = espressoTask.monitoringLog
    }
}
