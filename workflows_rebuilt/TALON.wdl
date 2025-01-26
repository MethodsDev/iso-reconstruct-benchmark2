version 1.0

task talonTask {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF

        Int cpu = 8
        Int numThreads = 16
        Int memoryGB = 128
        Int diskSizeGB = 1024
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/talon"

    }

    command <<<

        TALON-runner.py --output_prefix ~{sample_id} \
                        --genome ~{referenceGenomeFasta} \
                        --gtf ~{referenceAnnotationGTF}  \
                        --bam ~{inputBAM}
        
    >>>

    output {
        File talon_gtf = "~{sample_id}_talon.gtf"
        File talon_quant = "~{sample_id}_Quant_talon_abundance_filtered.tsv"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}

workflow talonWorkflow {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF
    }

    call talonTask {
        input:
            sample_id = sample_id,
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenomeFasta = referenceGenomeFasta,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotationGTF = referenceAnnotationGTF
    }

    output {
        File talon_gtf = talonTask.talon_gtf
        File talon_quant = talonTask.talon_quant
    }
}
