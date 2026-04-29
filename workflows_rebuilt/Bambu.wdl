version 1.0

# This task uses Bambu version 3.12.1
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
        String bambu_version = "3.12.1"
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/bambu:3.12.1"
        
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
        File bambu_quant = if (quant_only) then "~{sample_id}.quant_only.bambu-v~{bambu_version}.bambu.counts.txt" else if (defined(referenceAnnotationGTF)) then "~{sample_id}.Ref-Guided.bambu-v~{bambu_version}.bambu.counts.txt" else "~{sample_id}.Ref-Free.bambu-v~{bambu_version}.bambu.counts.txt"
        Array[File] bambu_gtf = glob("*.bambu.gtf")
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
            sample_id = sample_id,
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenomeFasta = referenceGenomeFasta,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotationGTF = referenceAnnotationGTF,
            quant_only = quant_only
    }

    output {
        File bambu_quant = bambuTask.bambu_quant
        Array[File] bambu_gtf = bambuTask.bambu_gtf
    }

}
