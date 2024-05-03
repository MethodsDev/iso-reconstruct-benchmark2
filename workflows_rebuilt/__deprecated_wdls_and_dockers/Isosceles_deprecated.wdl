version 1.0

task IsoscelesTask {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        String datasetName
        Int cpu = 16
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isosceles@sha256:f4e412b901dfb5a20da2c1aded7d66e2ca9e95e0291d20cd8bf21c70ec9fddb8"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
    }


    command <<<
        bash ~{monitoringScript} > monitoring.log &

        Rscript -<< "EOF"
        library(Isosceles)

        bam_file <- "~{inputBAM}"
        gtf_file <- "~{referenceAnnotation}"
        genome_fasta_file <- "~{referenceGenome}"
        
        bam_files <- c(Sample = bam_file)
        bam_parsed <- extract_read_structures(bam_files = bam_files)
        transcript_data <- prepare_transcripts(gtf_file = gtf_file, genome_fasta_file = genome_fasta_file, bam_parsed = bam_parsed, min_bam_splice_read_count = 2, min_bam_splice_fraction = 0.01)
        se_tcc <- prepare_tcc_se(bam_files = bam_files, transcript_data = transcript_data, run_mode = "de_novo_loose", min_read_count = 1, min_relative_expression = 0)
        se_transcript <- prepare_transcript_se(se_tcc = se_tcc, use_length_normalization = TRUE)
        export_gtf(se_transcript, "isoform_annotated.gtf")
        EOF

    >>>

    output {
        File isoscelesGTF = "isoform_annotated.gtf"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

workflow Isosceles {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        String datasetName
    }

    call IsoscelesTask {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName
    }

    output {
        File isoscelesGTF = IsoscelesTask.isoscelesGTF
        File monitoringLog = IsoscelesTask.monitoringLog
    }
}
