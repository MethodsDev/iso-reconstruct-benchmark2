version 1.0

# This task uses Salmon version 0.14.1
task salmonTask {
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
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/salmon@sha256:da4b4eb49e4e28e594e356938763569a1738ed0ea488f04e45fce2b3469c9db5"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
    }
    String OutDir = "Salmon_out"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        if [[ "~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            if [[ -z "~{referenceAnnotation_full}" ]]; then
                echo "Full annotation is not provided. Please provide the full annotation."
                exit 1
            fi  
            mkdir -p ~{OutDir}

            # Step 1: Convert BAM to FASTQ
            samtools bam2fq ~{inputBAM} > ~{OutDir}/tmp.fastq

            # Step 2: Create transcriptome from reference_full and reference genome using gffread
            gffread ~{referenceAnnotation_full} -g ~{referenceGenome} -w  ~{OutDir}/transcriptome.fa 

            # Step 3: Align reads to reference using minimap2
            minimap2 -a -t ~{numThreads} ~{OutDir}/transcriptome.fa  ~{OutDir}/tmp.fastq | samtools view --threads ~{numThreads} -bS | samtools sort --threads ~{numThreads} > ~{OutDir}/mapped.sorted.bam
            samtools index ~{OutDir}/mapped.sorted.bam
            
            cd ~{OutDir}
            salmon quant --noErrorModel --noLengthCorrection -t transcriptome.fa -l A -a mapped.sorted.bam -o quant -p ~{numThreads}
            mv ~{OutDir}/quant/quant.sf ~{OutDir}/Salmon_quant.sf
        fi
    >>>

    output {
        File salmonCounts = "~{OutDir}/Salmon_quant.sf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}


version 1.0

workflow salmonWorkflow {
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

    call salmonTask {
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
        File salmonCounts = salmonTask.salmonCounts
        File monitoringLog = salmonTask.monitoringLog
    }
}