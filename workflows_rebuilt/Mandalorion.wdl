version 1.0

# This task uses Mandalorion version 4.5.0
task mandalorionTask {
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
        Int diskSizeGB = 2048
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/mandalorion@sha256:60b8781d8a91e02124fc1a53806b4fb43bb45e8c9d044810fae2cd9bf99d9d65"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }
    
    String OutDir = "Mandalorion_out"

    command <<<
        bash ~{monitoringScript} > monitoring.log &
        mkdir -p ~{OutDir}

        if [ "~{ID_or_Quant_or_Both}" = "ID" ] || [ "~{ID_or_Quant_or_Both}" = "Both" ]; then
            samtools bam2fq ~{inputBAM} > ~{OutDir}/samtools.bam2fq.fastq

            python3 /usr/local/src/Mandalorion/Mando.py \
            -G ~{referenceGenome} \
            -f ~{OutDir}/samtools.bam2fq.fastq \
            -p ~{OutDir}_reffree \
            -t ~{numThreads}
            
            if [ -n "~{referenceAnnotation_reduced}" ]; then
                python3 /usr/local/src/Mandalorion/Mando.py \
                -G ~{referenceGenome} \
                -g ~{referenceAnnotation_reduced} \
                -f ~{OutDir}/samtools.bam2fq.fastq \
                -p ~{OutDir} \
                -t ~{numThreads}
                
                if [ -f ~{OutDir}/Isoforms.filtered.clean.gtf ]; then
                    mv ~{OutDir}/Isoforms.filtered.clean.gtf ~{OutDir}/Mandalorion_reduced.gtf             
                fi
            fi
            mv ~{OutDir}_reffree/Isoforms.filtered.clean.gtf ~{OutDir}/Mandalorion.gtf 
            if [ -d ~{OutDir}_reffree ]; then
                rm -r ~{OutDir}_reffree        
            fi
        fi
    >>>

    output {
        File? mandalorionReducedGTF = "~{OutDir}/Mandalorion_reduced.gtf"
        File？mandalorionGTF = "~{OutDir}/Mandalorion.gtf"
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

workflow mandalorionWorkflow {
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

    call mandalorionTask {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both,
    }

    output {
        File？mandalorionGTF = mandalorionTask.mandalorionGTF
        File? mandalorionReducedGTF = mandalorionTask.mandalorionReducedGTF
        File monitoringLog = mandalorionTask.monitoringLog
    }
}
