version 1.0

# This task uses FLAMES-py version 0.1
task flamesTask {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF
        String data_type
        
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 250
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/flames-py"
    }
    

    command <<<

        set -ex

        FLAMES-runner.py --output_prefix ~{sample_id} \
                         --genome ~{referenceGenomeFasta} \
                         --bam ~{inputBAM} \
                         --gtf ~{referenceAnnotationGTF} \
                         --data_type ~{data_type}
                
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

workflow flamesWorkflow {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF
        String data_type
    }

    call flamesTask {
        input:
            sample_id = sample_id,
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenomeFasta = referenceGenomeFasta,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotationGTF = referenceAnnotationGTF,
            data_type = data_type
    }

    output {
        File flames_gff3 = flamesTask.flames_gff3
        File flames_counts = flamesTask.flames_counts
    }

}
