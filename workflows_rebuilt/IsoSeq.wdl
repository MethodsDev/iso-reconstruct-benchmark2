version 1.0

# This task uses IsoSeq version 4.0.0 and Pigeon version 1.2.0
task isoseqTask {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File? referenceAnnotationGTF

        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 250
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoseq"
    }

    
    command <<<
        
        set -ex

        IsoSeq-runner.py --output_prefix ~{sample_id} \
                         --bam ~{inputBAM} \
                         --genome ~{referenceGenomeFasta} \
                         ~{"--gtf " + referenceAnnotationGTF} \
                         --ncpu ~{numThreads}

        
        
    >>>

    output {
        File isoseq_ref_free_GTF = "~{sample_id}.IsoSeq.ref-free.ID.gff"
        File? isoseq_ref_filtered_GTF = "~{sample_id}.IsoSeq.ref-filtered.ID.gff"
        File isoseq_counts = "~{sample_id}.IsoSeq.ref-free.ID.abundance.txt"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}

workflow isoseqWorkflow {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File? referenceAnnotationGTF
        
    }

    call isoseqTask {
        input:
            sample_id = sample_id,
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenomeFasta = referenceGenomeFasta,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotationGTF = referenceAnnotationGTF,
        
    }

    output {
        File isoseq_ref_free_GTF = isoseqTask.isoseq_ref_free_GTF
        File? isoseq_ref_filtered_GTF = isoseqTask.isoseq_ref_filtered_GTF
        File isoseq_counts = isoseqTask.isoseq_counts
    }
}
