version 1.0

task FilterTranscripts {
    input {
        String referenceGenome
        String gtf_path
        String expr_file_path
        String output_gtf_path
        Float threshold
        Int memoryGB
        Int diskSizeGB
        String docker
    }

    command <<<
        python /usr/local/bin/FilterTranscripts ~{gtf_path} ~{expr_file_path} ~{referenceGenome} ~{output_gtf_path} ~{threshold}
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 30
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }

    output {
        File filtered_gtf = output_gtf_path
    }
}

workflow TranscriptFiltering {
    input {
        String referenceGenome
        String gtf_path
        String expr_file_path
        String output_gtf_path
        Float threshold = 1.0
        Int memoryGB = 32
        Int diskSizeGB = 1024
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/filtertranscripts:latest"
    }

    call FilterTranscripts {
        input:
            referenceGenome = referenceGenome,
            gtf_path = gtf_path,
            expr_file_path = expr_file_path,
            output_gtf_path = output_gtf_path,
            threshold = threshold,
            docker = docker,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }

    output {
        File filtered_gtf = FilterTranscripts.filtered_gtf
    }
}
