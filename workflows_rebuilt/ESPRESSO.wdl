version 1.0

# This task uses ESPRESSO version 1.5.0
task espressoTask {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF

        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 256
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/espresso:latest"
    }


    command <<<

        set -ex

        ESPRESSO-runner.py --output_prefix ~{sample_id} \
                           --genome ~{referenceGenomeFasta} \
                           --gtf ~{referenceAnnotationGTF} \
                           --bam ~{inputBAM}

    >>>

    output {
        File espresso_gtf = "~{sample_id}.espresso.gtf"
        File espresso_counts = "~{sample_id}.espresso.counts.tsv"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}

workflow espressoWorkflow {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF
    }

    call espressoTask {
        input:
            sample_id = sample_id,
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenomeFasta = referenceGenomeFasta,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotationGTF = referenceAnnotationGTF
    }

    output {
        File espresso_gtf = espressoTask.espresso_gtf
        File espresso_counts = espressoTask.espresso_counts
    }
}
