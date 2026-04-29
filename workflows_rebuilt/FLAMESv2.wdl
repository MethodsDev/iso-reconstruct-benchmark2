version 1.0

# This task uses FLAMES v2 (mritchielab/FLAMES Bioconductor package)
task flamesV2Task {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF
        String data_type
        Boolean quant_only = false

        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 250
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/flames-v2:latest"
    }

    String quant_only_flag = if (quant_only) then "--quant_only" else ""

    command <<<
        set -ex

        FLAMESv2-runner.Rscript \
            --output_prefix ~{sample_id} \
            --genome ~{referenceGenomeFasta} \
            --bam ~{inputBAM} \
            --gtf ~{referenceAnnotationGTF} \
            --data_type ~{data_type} \
            --ncpu ~{numThreads} \
            ~{quant_only_flag}
    >>>

    output {
        File flames_gff3 = "~{sample_id}.FLAMES.gff3"
        File flames_counts = "~{sample_id}.FLAMES.counts.tsv"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}

workflow flamesV2Workflow {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF
        String data_type
        Boolean quant_only = false
    }

    call flamesV2Task {
        input:
            sample_id = sample_id,
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenomeFasta = referenceGenomeFasta,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotationGTF = referenceAnnotationGTF,
            data_type = data_type,
            quant_only = quant_only
    }

    output {
        File flames_gff3 = flamesV2Task.flames_gff3
        File flames_counts = flamesV2Task.flames_counts
    }
}
