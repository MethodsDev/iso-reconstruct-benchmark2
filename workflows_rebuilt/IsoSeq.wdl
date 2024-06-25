version 1.0

# This task uses IsoSeq version 4.0.0 and Pigeon version 1.2.0
task isoseqTask {
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
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoseq@sha256:c0d801579938050b15ee6e62adc31df8763b1d2aed0e4ff64e8026f154ab120e"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String OutDir = "IsoSeq_out"

    command <<<
        
        bash ~{monitoringScript} > monitoring.log &
        
        rm -rf ~{OutDir} && mkdir ~{OutDir}
    
        samtools bam2fq ~{inputBAM} > ~{OutDir}/temp.fastq
    
        mkdir ~{OutDir}/isoseq_reffree
        pbmm2 align --num-threads ~{numThreads} --preset ISOSEQ --sort ~{referenceGenome} ~{OutDir}/temp.fastq ~{OutDir}/isoseq_reffree/pbmm_aligned.bam
        isoseq3 collapse --do-not-collapse-extra-5exons ~{OutDir}/isoseq_reffree/pbmm_aligned.bam ~{OutDir}/isoseq_reffree/pbmm_aligned.gff
        cp ~{OutDir}/isoseq_reffree/pbmm_aligned.gff ~{OutDir}/IsoSeq.gff
    
        if [ -f "~{referenceAnnotation_reduced}" ]; then
            mkdir ~{OutDir}/isoseq
            pigeon prepare ~{OutDir}/isoseq_reffree/pbmm_aligned.gff
    
            cp ~{referenceAnnotation_reduced} ~{OutDir}/isoseq/ref_reduced.gtf
            cp ~{referenceGenome} ~{OutDir}/isoseq/ref_genome.fa
    
            pigeon prepare ~{OutDir}/isoseq/ref_reduced.gtf ~{OutDir}/isoseq/ref_genome.fa
    
            pigeon classify ~{OutDir}/isoseq_reffree/pbmm_aligned.sorted.gff ~{OutDir}/isoseq/ref_reduced.sorted.gtf ~{referenceGenome} --fl ~{OutDir}/isoseq_reffree/pbmm_aligned.flnc_count.txt -d ~{OutDir}/isoseq
            cp ~{OutDir}/isoseq_reffree/pbmm_aligned.sorted.gff ~{OutDir}/isoseq/pbmm_aligned.sorted.gff
            pigeon filter ~{OutDir}/isoseq/pbmm_aligned_classification.txt --isoforms ~{OutDir}/isoseq/pbmm_aligned.sorted.gff
            cp ~{OutDir}/isoseq/pbmm_aligned.sorted.filtered_lite.gff ~{OutDir}/IsoSeq_reduced.gff
        fi
    >>>

    output {
        File isoseqGTF = "~{OutDir}/IsoSeq.gff"
        File ?isoseqReducedGTF = "~{OutDir}/IsoSeq_reduced.gff"
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

workflow isoseqWorkflow {
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

    call isoseqTask {
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
        File isoseqGTF = isoseqTask.isoseqGTF
        File? isoseqReducedGTF = isoseqTask.isoseqReducedGTF
        File monitoringLog = isoseqTask.monitoringLog
    }
}
