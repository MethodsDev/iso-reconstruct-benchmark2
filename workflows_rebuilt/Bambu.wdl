version 1.0

# This task uses Bambu version 3.4.0
task bambuTask {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File? referenceAnnotationGTF
        Boolean quant_only
        
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 250
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/bambu:latest"
        
    }

    String quant_only_flag = if (quant_only) then "--quant_only" else ""

    command <<<

        set -ex

        bambu-runner.Rscript --bam ~{inputBAM} \
        --genome ~{referenceGenomeFasta} \
        ~{"--gtf " + referenceAnnotationGTF} \
        --output_prefix ~{sample_id} \
        ~{quant_only_flag}
            
    >>>

    output {
        File bambu_quant = if (quant_only) then "~{sample_id}.quant_only.bambu.counts.txt" else if (defined(referenceAnnotationGTF)) then "~{sample_id}.Ref-Guided.bambu.counts.txt" else "~{sample_id}.Ref-Free.bambu.counts.txt"
        File? bambu_gtf = if (defined(referenceAnnotationGTF)) then "~{sample_id}.Ref-Guided.bambu.gtf" else "~{sample_id}.Ref-Free.bambu.gtf"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}


workflow bambuWorkflow {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File? referenceAnnotationGTF
        Boolean quant_only
    }

    call bambuTask {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenomeFasta = referenceGenomeFasta,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotationGTF = referenceAnnotationGTF,
            quant_only = quant_only
    }

    output {
        File bambu_quant = bambuTask.bambu_quant
        File? bambu_gtf = bambuTask.bambu_gtf
    }

}

