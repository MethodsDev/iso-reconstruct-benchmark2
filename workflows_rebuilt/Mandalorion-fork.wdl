version 1.0

# This task uses Mandalorion version 4.5.0 (MDL fork version to be able to directly run on aligned bam files)
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
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 2048
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/mandalorion-fork@sha256:650bd77d9ec99064ad00eca5796d67fdd4aa4257bfccd996c88a749ee4d64322"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }
    
    String OutDir = "mandalorionfork_out"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        mkdir -p ~{OutDir}
        if [ "~{ID_or_Quant_or_Both}" = "ID" ] || [ "~{ID_or_Quant_or_Both}" = "Both" ]; then
            samtools bam2fq ~{inputBAM} > ~{OutDir}/samtools.bam2fq.fastq
            if samtools view -f 0x100 ~{inputBAM} | read -r; then
                samtools view -b -F 0x904 ~{inputBAM} > ~{OutDir}/output.bam
                samtools view -h -o samtools.view.sam ~{OutDir}/output.bam
            else
                samtools view -h -o samtools.view.sam ~{inputBAM}
            fi

#            python3 /usr/local/src/Mandalorion/Mando.py \
#            -G ~{referenceGenome} \
#            -f ~{OutDir}/samtools.bam2fq.fastq \
#            -p ~{OutDir}_reffree \
#            -t ~{numThreads} \
#            -s samtools.view.sam

            
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
        File? mandalorionforkGTF = "~{OutDir}/Mandalorionfork.gtf"
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
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }

    output {
        File? mandalorionforkGTF = mandalorionforkTask.mandalorionforkGTF
        File? mandalorionforkReducedGTF = mandalorionforkTask.mandalorionforkReducedGTF
        File monitoringLog = mandalorionforkTask.monitoringLog
    }
}
