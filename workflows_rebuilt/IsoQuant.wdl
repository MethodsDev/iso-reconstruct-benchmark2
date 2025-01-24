version 1.0

task isoquantTask {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File? referenceAnnotationGTF
        Boolean quant_only
        String data_type
        
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:latest"
    }

    String quant_only_flag = if (quant_only) then "--quant_only" else ""
    
    
    command <<<

        set -ex

        IsoQuant-runner.py  --genome ~{referenceGenomeFasta} \
                            --bam ~{inputBAM} \
                            --gtf ~{referenceAnnotationGTF} \
                            --data_type ~{data_type} \
                            --output_prefix ~{sample_id} \
                            ~{quant_only_flag}


        tar -czhf ~{sample_id}.IsoQuant_outdir.tar.gz isoquant_output_dir/ 

        
    >>>

    output {
        File? isoquant_gtf = "~{sample_id}.IsoQuant.gtf"
        File isoquant_counts = "~{sample_id}.IsoQuant.counts.tsv"
        File isoquant_outdir_tgz = "~{sample_id}.IsoQuant_outdir.tar.gz"

    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}

workflow isoquantWorkflow {
    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex
        File? referenceAnnotationGTF
        Boolean quant_only
        String data_type
    }

    call isoquantTask {
        input:
            sample_id = sample_id,
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenomeFasta = referenceGenomeFasta,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotationGTF = referenceAnnotationGTF,
            quant_only = quant_only,
            data_type = data_type
    }

    output {
        File? isoquant_gtf = isoquantTask.isoquant_gtf
        File isoquant_counts = isoquantTask.isoquant_counts
        File isoquant_outdir_tgz = isoquantTask.isoquant_outdir_tgz
    }
}
