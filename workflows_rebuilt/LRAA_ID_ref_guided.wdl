version 1.0

task splitBAMByChromosome {
    input {
        File inputBAM
        String main_chromosomes
        String docker
        Int threads
        File referenceGenome
        File referenceAnnotation_reduced
        Int memoryGB
        Int diskSizeGB
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        set -eo pipefail
        mkdir -p split_bams
        
        # Check if BAM index exists, if not, index the input BAM
        if [ ! -f "~{inputBAM}.bai" ]; then
            samtools index -@ ~{threads} ~{inputBAM}
        fi
        
        # Loop through each chromosome
        for chr in ~{main_chromosomes}; do
            # Generate chromosome-specific BAM
            samtools view -@ ~{threads} -b ~{inputBAM} $chr > split_bams/$chr.bam
            
            # Generate chromosome-specific FASTA from the whole genome
            samtools faidx ~{referenceGenome} $chr > split_bams/$chr.genome.fasta
            
            # Generate chromosome-specific GTF for reduced annotation, if available
            if [ -f "~{referenceAnnotation_reduced}" ]; then
                cat ~{referenceAnnotation_reduced} | perl -lane 'if ($F[0] eq "'$chr'") { print; }' > split_bams/$chr.reduced.annot.gtf
            fi
        done
    >>>

    output {
        Array[File] chromosomeBAMs = glob("split_bams/*.bam")
        Array[File] chromosomeFASTAs = glob("split_bams/*.genome.fasta")
        Array[File] reducedAnnotations = glob("split_bams/*.reduced.annot.gtf")
    }

    runtime {
        docker: docker
        bootDiskSizeGb: 30
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task lraaPerChromosome {
    input {
        File inputBAM
        File referenceGenome
        String OutDir
        String docker
        Int numThreads
        Boolean? LRAA_no_norm
        File referenceAnnotation_reduced
        Int memoryGB
        Int diskSizeGB
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String no_norm_flag = if defined(LRAA_no_norm) && LRAA_no_norm then "--no_norm" else ""
    
    command <<<
        bash ~{monitoringScript} > monitoring.log &

        mkdir -p ~{OutDir}/ID_reduced
    
        /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/ID_reduced/LRAA_reduced \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_reduced} --CPU 1
    >>>
    
    output {
        File lraaID_reduced_GTF = "~{OutDir}/ID_reduced/LRAA_reduced.gtf"
    }
    runtime {
        docker: docker
        bootDiskSizeGb: 30
        cpu: "~{numThreads}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task mergeResults {
    input {
        Array[File] reducedGtfFiles
        String outputFilePrefix
        String docker
        Int memoryGB
        Int diskSizeGB
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        set -eo pipefail
    
        # Initialize output file
        reduced_gtf_output="~{outputFilePrefix}_merged_reduced.gtf"
        
        # Directly concatenate all input files into the output file
        cat ~{sep=' ' reducedGtfFiles} > "$reduced_gtf_output"
    >>>

    output {
        File mergedReducedGtfFile = "~{outputFilePrefix}_merged_reduced.gtf"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

workflow lraaWorkflow {
    input {
        File? inputBAM
        Array[File]? inputBAMArray
        File referenceGenome
        Int numThreads = 4
        Int memoryGB = 32
        Int diskSizeGB = 1024
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        File referenceAnnotation_reduced
        String main_chromosomes = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
        Boolean? LRAA_no_norm
    }

    String OutDir = "LRAA_out"

    Array[File] chromosomeBAMs
    Array[File] chromosomeFASTAs
    Array[File] reducedAnnotations

    if (defined(inputBAM)) {
        call splitBAMByChromosome {
            input:
                inputBAM = inputBAM,
                main_chromosomes = main_chromosomes,
                docker = docker,
                threads = numThreads,
                referenceGenome = referenceGenome,
                referenceAnnotation_reduced = referenceAnnotation_reduced,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
        }

        chromosomeBAMs = splitBAMByChromosome.chromosomeBAMs
        chromosomeFASTAs = splitBAMByChromosome.chromosomeFASTAs
        reducedAnnotations = splitBAMByChromosome.reducedAnnotations
    } else {
        chromosomeBAMs = inputBAMArray
        chromosomeFASTAs = [referenceGenome] # Assuming the reference genome is the same for all chromosomes
        reducedAnnotations = [referenceAnnotation_reduced] # Assuming the reduced annotation is the same for all chromosomes
    }

    scatter (i in range(length(chromosomeBAMs))) {
        call lraaPerChromosome {
            input:
                inputBAM = chromosomeBAMs[i],
                referenceGenome = chromosomeFASTAs[i],
                OutDir = OutDir,
                docker = docker,
                numThreads = numThreads,
                LRAA_no_norm = LRAA_no_norm,
                referenceAnnotation_reduced = reducedAnnotations[i],
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
        }
    }

    call mergeResults {
        input:
            reducedGtfFiles = lraaPerChromosome.lraaID_reduced_GTF,
            outputFilePrefix = "merged",
            docker = docker,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }
    
    output {
        File mergedReducedGTF = mergeResults.mergedReducedGtfFile
        Array[File]? splitBAMs = if defined(inputBAM) then splitBAMByChromosome.chromosomeBAMs else []
    }
}
