version 1.0

task isoscelesTask {
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
        Int diskSizeGB = 250
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoceles2:latest"
    }

    String quant_only_str = if (quant_only) then "--quant_only" else ""
    
    command <<<
        set -ex
        
        Isosceles-runner.Rscript \
            --bam ~{inputBAM} \
            --genome ~{referenceGenomeFasta} \
            --ncpu ~{numThreads} \
            --gtf ~{referenceAnnotationGTF} \
            ~{quant_only_str} \
             --output_prefix ~{sample_id}
                
    >>>
    
    output {
        File isosceles_counts = if (quant_only) then "~{sample_id}.isosceles.counts" else "~{sample_id}.isosceles.refguided_novel_ID.counts"
        File? isosceles_gtf = "~{sample_id}.isosceles.refguided_novel_ID.gtf"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}

workflow isoscelesWorkflow {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF
        Boolean quant_only
    }

    call isoscelesTask {
        input:
            sample_id = sample_id,
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenomeFasta = referenceGenomeFasta,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotationGTF = referenceAnnotationGTF,
            quant_only = quant_only
    }

    output {
        File isosceles_counts = isoscelesTask.isosceles_counts
        File? isosceles_gtf = isoscelesTask.isosceles_gtf
    }
}
