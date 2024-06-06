version 1.0

# This task uses LRAA version 0.0.8
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
        Int? min_mapping_quality
        Boolean? LRAA_no_norm
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String OutDir = "LRAA_out"
    String min_mapping_quality_flag = if (defined(min_mapping_quality)) then "--min_mapping_quality=" + min_mapping_quality else ""

    String no_norm_flag = if (defined(LRAA_no_norm) && LRAA_no_norm) then "--no_norm" else ""


    command <<<
        bash ~{monitoringScript} > monitoring.log &
        
        mkdir -p ~{OutDir}/ID ~{OutDir}/ID_reduced ~{OutDir}/Quant ~{OutDir}/Quant_noEM

        if [[ "~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/ID/LRAA \
                                 ~{no_norm_flag}

        fi

        if [[ ("~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both") && -n "~{referenceAnnotation_reduced}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/ID_reduced/LRAA_reduced \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_reduced}


        fi

        if [[ ("~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both") && -n "~{referenceAnnotation_full}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/Quant/LRAA \
                                 --quant_only \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_full} \
                                 --EM


            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/Quant_noEM/LRAA.noEM \
                                 --quant_only \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_full}
        fi

        if [[ ("~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both") && -n "~{referenceAnnotation_full}" && -n "~{min_mapping_quality}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/Quant_noEM_minMapQ/LRAA.noEM.minMapQ \
                                 --quant_only \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_full} \
                                 ~{min_mapping_quality_flag}
        fi
    >>>

    output {
        File? lraaGTF = "~{OutDir}/ID/LRAA.gtf"
        File? lraaReducedGTF = "~{OutDir}/ID_reduced/LRAA_reduced.gtf"
        File? lraaCounts = "~{OutDir}/Quant/LRAA.quant.expr"
        File? lraaCounts_noEM = "~{OutDir}/Quant_noEM/LRAA.noEM.quant.expr"
        File? lraa_quant_tracking = "~{OutDir}/Quant/LRAA.quant.tracking"
        File? lraa_quant_tracking_noEM = "~{OutDir}/Quant_noEM/LRAA.noEM.quant.tracking"
        File? lraaCounts_noEM_minMapQ = "~{OutDir}/Quant_noEM_minMapQ/LRAA.noEM.minMapQ.quant.expr"
        File? lraa_quant_tracking_noEM_minMapQ = "~{OutDir}/Quant_noEM_minMapQ/LRAA.noEM.minMapQ.quant.tracking"
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
        Int? min_mapping_quality
        Boolean? LRAA_no_norm
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
            ID_or_Quant_or_Both = ID_or_Quant_or_Both,
            min_mapping_quality = min_mapping_quality,
            LRAA_no_norm = LRAA_no_norm
    }

    output {
        File? lraaGTF = lraaTask.lraaGTF
        File? lraaReducedGTF = lraaTask.lraaReducedGTF
        File? lraaCounts = lraaTask.lraaCounts
        File? lraaCounts_noEM = lraaTask.lraaCounts_noEM
        File? lraa_quant_tracking = lraaTask.lraa_quant_tracking
        File? lraa_quant_tracking_noEM = lraaTask.lraa_quant_tracking_noEM
        File? lraaCounts_noEM_minMapQ = lraaTask.lraaCounts_noEM_minMapQ
        File? lraa_quant_tracking_noEM_minMapQ = lraaTask.lraa_quant_tracking_noEM_minMapQ
        File monitoringLog = lraaTask.monitoringLog
    }
}
