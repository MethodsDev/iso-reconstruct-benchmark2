version 1.0

# This task uses FLAMES-py version 0.1
task flamesTask {
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
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/flames-py@sha256:3ce1f2c8c20088945ce885dcbf486ae377873e30975dad55e54235c3873648d4"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
    }
    
    String OutDir = "FLAMES_out"

    command <<<
        bash ~{monitoringScript} > monitoring.log &
        mkdir -p ~{OutDir}

        if [ "~{dataType}" = "pacbio_ccs" ]; then
            cat > ~{OutDir}/config.json << EOF
            {
                "comment":"this is the default config for nanopore single cell long read data using 10X RNA-seq kit. use splice annotation in alignment.",
                "pipeline_parameters":{
                    "do_genome_alignment":true,
                    "do_isoform_identification":true,
                    "do_read_realignment":true,
                    "do_transcript_quantification":true
                },
                "global_parameters":{
                    "generate_raw_isoform":false,
                    "has_UMI":true
                },
                "isoform_parameters":{
                    "MAX_DIST":10,
                    "MAX_TS_DIST":120,
                    "MAX_SPLICE_MATCH_DIST":10,
                    "min_fl_exon_len":40,
                    "Max_site_per_splice":3,
                    "Min_sup_cnt":5,
                    "Min_cnt_pct":0.001,
                    "Min_sup_pct":0.2,
                    "strand_specific":1,
                    "remove_incomp_reads":4,
                    "random_seed":666666
                },
                "alignment_parameters":{
                    "use_junctions":true,
                    "no_flank":false,
                    "seed":2022
                },
                "realign_parameters":{
                    "use_annotation":true
                },
                "transcript_counting":{
                    "min_tr_coverage":0.4,
                    "min_read_coverage":0.4
                }
            }
EOF
        else
            cat > ~{OutDir}/config.json << EOF
            {
                "comment":"this is the default config for nanopore single cell long read data using 10X RNA-seq kit. use splice annotation in alignment.",
                "pipeline_parameters":{
                    "do_genome_alignment":true,
                    "do_isoform_identification":true,
                    "do_read_realignment":true,
                    "do_transcript_quantification":true
                },
                "global_parameters":{
                    "generate_raw_isoform":false,
                    "has_UMI":true
                },
                "isoform_parameters":{
                    "MAX_DIST":10,
                    "MAX_TS_DIST":120,
                    "MAX_SPLICE_MATCH_DIST":10,
                    "min_fl_exon_len":40,
                    "Max_site_per_splice":3,
                    "Min_sup_cnt":5,
                    "Min_cnt_pct":0.001,
                    "Min_sup_pct":0.2,
                    "strand_specific":0,
                    "remove_incomp_reads":4,
                    "random_seed":666666
                },
                "alignment_parameters":{
                    "use_junctions":true,
                    "no_flank":false,
                    "seed":2022
                },
                "realign_parameters":{
                    "use_annotation":true
                },
                "transcript_counting":{
                    "min_tr_coverage":0.4,
                    "min_read_coverage":0.4
                }
            }
EOF
        fi

        if [ "~{ID_or_Quant_or_Both}" = "ID" -o "~{ID_or_Quant_or_Both}" = "Both" ] && [ -n "~{referenceAnnotation_reduced}" ]; then
            rm -rf ~{OutDir}/fq && mkdir ~{OutDir}/fq
            samtools bam2fq ~{inputBAM} > ~{OutDir}/fq/temp.fastq

            python3 /usr/local/src/FLAMES/python/bulk_long_pipeline.py \
            --gff3 ~{referenceAnnotation_reduced} \
            --genomefa ~{referenceGenome} \
            --fq_dir ~{OutDir}/fq \
            --inbam ~{inputBAM} \
            --outdir ~{OutDir} \
            --config_file ~{OutDir}/config.json
            
            mv ~{OutDir}/isoform_annotated.gff3 ~{OutDir}/FLAMES_reduced.gff3
        fi
    >>>

    output {
        File? flamesReducedGTF = "~{OutDir}/FLAMES_reduced.gff3"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}

workflow flamesWorkflow {
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

    call flamesTask {
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
        File? flamesReducedGTF = flamesTask.flamesReducedGTF
        File monitoringLog = flamesTask.monitoringLog
    }
}
