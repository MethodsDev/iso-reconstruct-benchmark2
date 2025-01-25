version 1.0

# This task uses Oarfish version 4.5.0
task oarfishTask {
    input {
        String sample_id
        File inputBAM # will convert to fastq unless inputFASTQ provided 
        File? inputFASTQ # if fastq provided, uses it instead of the bam
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF
        String oarfish_seq_tech # ont-cdna, ont-drna, pac-bio, or pac-bio-hifi
        String oarfish_mode # byAlignment or byReads
        
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 250
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/oarfish"

    }


    command <<<

        set -ex

        inputFASTQ=~{inputFASTQ}
        if [[ "$inputFASTQ" == "" ]]; then
           samtools bam2fq ~{inputBAM} > reads.fastq
           inputFASTQ="reads.fastq"
        fi
        
        Oarfish_runner.py \
                     --output_prefix ~{sample_id} \
                     --genome ~{referenceGenomeFasta} \
                     --gtf ~{referenceAnnotationGTF} \
                     --fastq $inputFASTQ \
                     --threads ~{cpu} \
                     --mode ~{oarfish_mode} \
                     --seq_tech ~{oarfish_seq_tech}
        
    >>>

    output {
        File oarfish_quant = "~{sample_id}.Oarfish.~{oarfish_mode}.quant"
       
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}


workflow oarfishWorkflow {
    input {
        String sample_id
        File inputBAM # will convert to fastq unless inputFASTQ provided
        File? inputFASTQ  # if fastq provided, uses it instead of the bam
        File referenceGenomeFasta
        File referenceGenomeIndex
        File referenceAnnotationGTF
        String oarfish_seq_tech # ont-cdna, ont-drna, pac-bio, or pac-bio-hifi
        String oarfish_mode # byAlignment or byReads 
    }

    call oarfishTask {
        input:
            sample_id = sample_id,
            inputBAM = inputBAM,
            inputFASTQ = inputFASTQ,
            referenceGenomeFasta = referenceGenomeFasta,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotationGTF = referenceAnnotationGTF,
            oarfish_seq_tech = oarfish_seq_tech,
            oarfish_mode = oarfish_mode
    }

    output {
        File oarfish_quant = oarfishTask.oarfish_quant
    }
}
