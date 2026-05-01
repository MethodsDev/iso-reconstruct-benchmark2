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
        String isoquant_version_tag = "v3.13.0"
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:latest"
    }

    String quant_only_flag = if (quant_only) then "--quant_only" else ""
    
    
    command <<<

        set -ex

        RUNNER_STATUS=0
        IsoQuant-runner.py  --genome ~{referenceGenomeFasta} \
                            --bam ~{inputBAM} \
                            ~{"--gtf " + referenceAnnotationGTF} \
                            --data_type ~{data_type} \
                            --output_prefix ~{sample_id} \
                            --isoquant_version_tag ~{isoquant_version_tag} \
                            ~{quant_only_flag} || RUNNER_STATUS=$?

        OUTDIR="~{sample_id}.isoquant.outdir/OUT"
        COUNTS_DEST="~{sample_id}.isoquant-~{isoquant_version_tag}.IsoQuant.counts.tsv"
        GTF_DEST="~{sample_id}.isoquant-~{isoquant_version_tag}.IsoQuant.gtf"

        # Fallback for when the runner inside the docker image fails to
        # produce COUNTS_DEST (e.g., because a stale runner only knows
        # the old `transcript_model_counts.tsv` candidate and that file
        # is missing in newer IsoQuant versions). Order matters: in
        # ref-guided mode v3.13.0 emits both `transcript_counts.tsv`
        # (reference transcripts only) and
        # `discovered_transcript_counts.tsv` (reference + newly
        # discovered models). Always prefer the discovered file so the
        # novels make it through to downstream benchmarking.
        if [ ! -s "${COUNTS_DEST}" ]; then
            if [ -s "${OUTDIR}/OUT.transcript_model_counts.tsv" ]; then
                cp "${OUTDIR}/OUT.transcript_model_counts.tsv" "${COUNTS_DEST}"
            elif [ -s "${OUTDIR}/OUT.discovered_transcript_counts.tsv" ]; then
                cp "${OUTDIR}/OUT.discovered_transcript_counts.tsv" "${COUNTS_DEST}"
            elif [ -s "${OUTDIR}/OUT.transcript_counts.tsv" ]; then
                cp "${OUTDIR}/OUT.transcript_counts.tsv" "${COUNTS_DEST}"
            fi
        fi

        if [ ! -s "${GTF_DEST}" ] && [ -s "${OUTDIR}/OUT.transcript_models.gtf" ]; then
            cp "${OUTDIR}/OUT.transcript_models.gtf" "${GTF_DEST}"
        fi

        # Only ignore runner failures when required outputs are successfully recovered.
        if [ "${RUNNER_STATUS}" -ne 0 ] && [ ! -s "${COUNTS_DEST}" ]; then
            exit "${RUNNER_STATUS}"
        fi

        if [ -d "isoquant_output_dir" ]; then
            tar -czhf ~{sample_id}.IsoQuant_outdir.tar.gz isoquant_output_dir/
        else
            tar -czhf ~{sample_id}.IsoQuant_outdir.tar.gz "~{sample_id}.isoquant.outdir/"
        fi

        
    >>>

    output {
        File? isoquant_gtf = "~{sample_id}.isoquant-~{isoquant_version_tag}.IsoQuant.gtf"
        File isoquant_counts = "~{sample_id}.isoquant-~{isoquant_version_tag}.IsoQuant.counts.tsv"
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
