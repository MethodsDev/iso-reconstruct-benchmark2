version 1.0

# This task uses LRAA version 0.0.4
task lraaTask {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String dataType
        String ID_or_Quant_or_Both
        Boolean? LRAA_no_norm = false
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String OutDir = "LRAA_out"
    Boolean quant_only = false
    Boolean? no_norm = false

    command <<<
        bash ~{monitoringScript} > monitoring.log &
        
        rm -rf ~{OutDir} && mkdir ~{OutDir} 

        out_prefix=lraa

        /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/~{ID_or_Quant_or_Both} \
                                 ~{true='--quant_only' false='' quant_only} \
                                 ~{true="--no_norm" false="" LRAA_no_norm} \
                                 ~{"--gtf " + referenceAnnotation_full} \
                                 --output_prefix ~{OutDir}/~{ID_or_Quant_or_Both}${out_prefix}
    >>>

    output {
        File lraaGTF = "~{OutDir}/~{ID_or_Quant_or_Both}.gtf"
        File? lraaReducedGTF = "~{OutDir}/~{ID_or_Quant_or_Both}_reduced.gtf"
        File? isoquantCounts = "~{OutDir}/~{ID_or_Quant_or_Both}_quant.tsv"
        File? isoquantGTF_with_polyA = "~{OutDir}/~{ID_or_Quant_or_Both}_with_polyA.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

workflow lraaWorkflow {
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

    call lraaTask {
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
        File lraaGTF = lraaTask.lraaGTF
        File? lraaReducedGTF = lraaTask.lraaReducedGTF
        File? isoquantCounts = lraaTask.isoquantCounts
        File? isoquantGTF_with_polyA = lraaTask.isoquantGTF_with_polyA
        File monitoringLog = lraaTask.monitoringLog
    }
}
