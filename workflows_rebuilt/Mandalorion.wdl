version 1.0

task mandalorionTask {
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
        Int diskSizeGB = 2048
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/mandalorion"

    }
    


    command <<<
        
        set -ex

        Mandalorian-runner.py --output_prefix  ~{sample_id} \
                              --genome ~{referenceGenomeFasta} \
                              --bam ~{inputBAM} \
                              ~{"--gtf " + referenceAnnotationGTF} \
                              --delayTime 60

    >>>

    output {
        File mandalorion_gtf = "~{sample_id}.Mandalorian.Isoforms.filtered.clean.gtf"
        File mandalorion_quant = "~{sample_id}.Mandalorian.Isoforms.filtered.clean.quant"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}

workflow mandalorionWorkflow {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File? referenceAnnotationGTF
        
    }

    call mandalorionTask {
        input:
            sample_id = sample_id,
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenomeFasta = referenceGenomeFasta,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotationGTF = referenceAnnotationGTF,
    }

    output {
        File mandalorion_gtf = mandalorionTask.mandalorion_gtf
        File mandalorion_quant = mandalorionTask.mandalorion_quant
    }
}

