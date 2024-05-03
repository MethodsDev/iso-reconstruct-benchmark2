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
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/espresso@sha256:f538303f6457c55e7b3c2a45081e6d8e3053e6f76e56bc65631b7f4aa290b026"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String OutDir = "ESPRESSO_out"
    String samples_filename = "espresso_samples.tsv"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        # Convert BAM to SAM
        samtools view -h -o input.sam ~{inputBAM}

        # Create espresso_samples.tsv
        echo -e "input.sam\tespresso" > ~{OutDir}/{samples_filename}

        if [[ "~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both" || "~{ID_or_Quant_or_Both}" == "ID" && -n "~{referenceAnnotation_reduced}" ]]; then
            mkdir -p ~{OutDir}
            perl /usr/src/app/espresso/src/ESPRESSO_S.pl --sort_buffer_size 16G -L ~{OutDir}/{samples_filename} -F ~{referenceGenome} -A ~{referenceAnnotation_full} -O ~{OutDir} -T ~{numThreads}
            perl /usr/src/app/espresso/src/ESPRESSO_C.pl --sort_buffer_size 16G -I ~{OutDir} -F ~{referenceGenome} -X 0 -T ~{numThreads}
            perl /usr/src/app/espresso/src/ESPRESSO_Q.pl -L ~{OutDir}/espresso_samples.tsv.updated -A ~{referenceAnnotation_full} -T ~{numThreads}
            mv ~{OutDir}/espresso_samples_N2_R0_abundance.esp ~{OutDir}/ESPRESSO_quant.txt
            if [[ "~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
                mv ~{OutDir}/espresso_samples_N2_R0_updated.gtf ~{OutDir}/ESPRESSO_reduced.gtf
            fi
        else
            echo "Reference annotation reduced is not provided. Cannot run in ID mode."
        fi
    >>>

    output {
        File? espressoCounts = "~{OutDir}/ESPRESSO_quant.txt"
        File? espressoReducedGTF = "~{OutDir}/ESPRESSO_reduced.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: 16
        memory: "256 GiB"
        disks: "local-disk 500 HDD"
        docker: docker
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
        Int cpu
        Int numThreads
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript
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
            ID_or_Quant_or_Both = ID_or_Quant_or_Both,
            cpu = cpu,
            numThreads = numThreads,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB,
            docker = docker,
            monitoringScript = monitoringScript
    }

    output {
        File? espressoCounts = espressoTask.espressoCounts
        File? espressoReducedGTF = espressoTask.espressoReducedGTF
        File monitoringLog = espressoTask.monitoringLog
    }
}
