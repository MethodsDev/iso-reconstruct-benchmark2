version 1.0

task stringtieTask {
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
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/stringtie"
    }

    String quant_only_flag = if (quant_only) then "--quant_only" else ""

    command <<<

        set -ex

        stringtie-runner.py --genome ~{referenceGenomeFasta} \
                            --bam ~{inputBAM} \
                            ~{"--gtf " + referenceAnnotationGTF} \
                            ~{quant_only_flag} \
                            --output_prefix ~{sample_id}

    >>>
    
    output {

        File stringtie_gtf = "~{sample_id}.stringtie.gtf"
        File stringtie_quant = "~{sample_id}.stringtie.quant.tsv"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}

workflow stringtieWorkflow {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File? referenceAnnotationGTF
        Boolean quant_only
    }

    call stringtieTask {
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
        File stringtie_gtf = stringtieTask.stringtie_gtf
        File stringtie_quant = stringtieTask.stringtie_quant
    }
}
