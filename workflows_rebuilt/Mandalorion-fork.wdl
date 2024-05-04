version 1.0

# This task uses Mandalorion version 4.5.0 (MDL fork version)
task mandalorionforkTask {
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
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/mandalorion-fork@sha256:60b8781d8a91e02124fc1a53806b4fb43bb45e8c9d044810fae2cd9bf99d9d65"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
    }
    
    String OutDir = "mandalorionfork_out"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        if [ "~{ID_or_Quant_or_Both}" = "ID" ] || [ "~{ID_or_Quant_or_Both}" = "Both" ]; then
            samtools bam2fq ~{inputBAM} > ~{OutDir}/samtools.bam2fq.fastq
            samtools view -F 0x904 ~{inputBAM} > ~{OutDir}/output.bam

            samtools view -h -o samtools.view.sam ~{OutDir}/output.bam

            python3 /usr/local/src/Mandalorion/Mando.py \
            -G ~{referenceGenome} \
            -f ~{OutDir}/samtools.bam2fq.fastq \
            -p ~{OutDir}_reffree \
            -t ~{numThreads} \
            -s samtools.view.sam

            
            if [ -n "~{referenceAnnotation_reduced}" ]; then
                python3 /usr/local/src/Mandalorion/Mando.py \
                -G ~{referenceGenome} \
                -g ~{referenceAnnotation_reduced} \
                -f ~{OutDir}/samtools.bam2fq.fastq \
                -p ~{OutDir} \
                -t ~{numThreads} \
                -s samtools.view.sam

                
                if [ -f ~{OutDir}/Isoforms.filtered.clean.gtf ]; then
                    mv ~{OutDir}/Isoforms.filtered.clean.gtf ~{OutDir}/Mandalorionfork_reduced.gtf             
                fi
            fi
            mv ~{OutDir}_reffree/Isoforms.filtered.clean.gtf ~{OutDir}/Mandalorionfork.gtf 
            if [ -d ~{OutDir}_reffree ]; then
                rm -r ~{OutDir}_reffree        
            fi
        fi
    >>>

    output {
        File? mandalorionforkReducedGTF = "~{OutDir}/Mandalorionfork_reduced.gtf"
        File mandalorionforkGTF = "~{OutDir}/Mandalorionfork.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

workflow mandalorionforkWorkflow {
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

    call mandalorionforkTask {
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
        File mandalorionforkGTF = mandalorionforkTask.mandalorionforkGTF
        File? mandalorionforkReducedGTF = mandalorionforkTask.mandalorionforkReducedGTF
        File monitoringLog = mandalorionforkTask.monitoringLog
    }
}
