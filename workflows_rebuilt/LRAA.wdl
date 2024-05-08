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
        Boolean? LRAA_no_norm
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String OutDir = "LRAA_out"
    String out_prefix = "LRAA"
    String no_norm_flag = if (defined(LRAA_no_norm) && LRAA_no_norm) then "--no_norm" else ""

    command <<<
        bash ~{monitoringScript} > monitoring.log &
        
        mkdir -p ~{OutDir}/ID ~{OutDir}/ID_reduced ~{OutDir}/Quant ~{OutDir}/Quant_noEM

        if [[ "~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            /home/jupyter/tools/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/ID/~{out_prefix} \
                                 ~{no_norm_flag}
        fi

        if [[ ("~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both") && -n "~{referenceAnnotation_reduced}" ]]; then
            /home/jupyter/tools/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/ID_reduced/~{out_prefix}_reduced \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_reduced}
        fi

        if [[ ("~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both") && -n "~{referenceAnnotation_full}" ]]; then
            if [[ -n "~{referenceAnnotation_reduced}" ]]; then
                /home/jupyter/tools/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/Quant/~{out_prefix} \
                                     --quant_only \
                                     ~{no_norm_flag} \
                                     --gtf ~{referenceAnnotation_full}

                /home/jupyter/tools/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/Quant_noEM/~{out_prefix}.noEM \
                                     --quant_only \
                                     ~{no_norm_flag} \
                                     --gtf ~{referenceAnnotation_full} \
                                     --no_EM
            fi
        fi
    >>>

    output {
        File lraaGTF = "~{OutDir}/ID/~{out_prefix}.gtf"
        File? lraaReducedGTF = "~{OutDir}/ID_reduced/~{out_prefix}_reduced.gtf"
        File? lraaCounts = "~{OutDir}/Quant/~{out_prefix}.quant.expr"
        File? lraaCounts_noEM = "~{OutDir}/Quant_noEM/~{out_prefix}.noEM.quant.expr"
        File? lraa_quant_tracking = "~{OutDir}/Quant/~{out_prefix}.quant.tracking"
        File? lraa_quant_tracking_noEM = "~{OutDir}/Quant_noEM/~{out_prefix}.noEM.quant.tracking"
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
            LRAA_no_norm = LRAA_no_norm,
    }

    output {
        File lraaGTF = lraaTask.lraaGTF
        File? lraaReducedGTF = lraaTask.lraaReducedGTF
        File? lraaCounts = lraaTask.lraaCounts
        File? lraaCounts_noEM = lraaTask.lraaCounts_noEM
        File? lraa_quant_tracking = lraaTask.lraa_quant_tracking
        File? lraa_quant_tracking_noEM = lraaTask.lraa_quant_tracking_noEM
        File monitoringLog = lraaTask.monitoringLog
    }
}
