version 1.0

# This task uses Flair version 3.0.0
task flairTask {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF
        Boolean quant_only
        
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 512
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/flair:latest"
        String flair_version_tag = "v3.0.0"
    }

    String quant_only_flag = if (quant_only) then "--quant_only" else ""

    command <<<
        
        set -ex

        flair-runner.py --genome ~{referenceGenomeFasta} \
                        --bam ~{inputBAM} \
                        --gtf ~{referenceAnnotationGTF} \
                        --output_prefix ~{sample_id} \
                        --flair_version_tag ~{flair_version_tag} \
                        ~{quant_only_flag}
        
    >>>

    output {
        File? flair_gtf = "~{sample_id}.flair-~{flair_version_tag}.flair.gtf"
        File flair_counts = "~{sample_id}.flair-~{flair_version_tag}.counts.tsv"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}


workflow flairWorkflow {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF
        Boolean quant_only
        String flair_version_tag = "v3.0.0"
    }

    call flairTask {
        input:
            sample_id = sample_id,
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenomeFasta = referenceGenomeFasta,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotationGTF = referenceAnnotationGTF,
            quant_only = quant_only,
            flair_version_tag = flair_version_tag
    }

    output {
        File? flair_gtf = flairTask.flair_gtf
        File flair_counts = flairTask.flair_counts
    }
}
