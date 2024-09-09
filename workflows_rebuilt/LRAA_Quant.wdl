version 1.0

task splitBAMByChromosome {
    input {
        File inputBAM
        String main_chromosomes
        String docker
        Int threads
        File referenceGenome
        File referenceAnnotation_full
        Int memoryGB
        Int diskSizeGB
    }

    command <<<
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
            
            # Generate chromosome-specific GTF for full annotation, if available
            if [ -f "~{referenceAnnotation_full}" ]; then
                cat ~{referenceAnnotation_full} | perl -lane 'if ($F[0] eq "'$chr'") { print; }' > split_bams/$chr.full.annot.gtf
            fi
        done
    >>>

    output {
        Array[File] chromosomeBAMs = glob("split_bams/*.bam")
        Array[File] chromosomeFASTAs = glob("split_bams/*.genome.fasta")
        Array[File] fullAnnotations = glob("split_bams/*.full.annot.gtf")
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
        Int? LRAA_min_mapping_quality
        File referenceAnnotation_full
        Int memoryGB
        Int diskSizeGB
    }

    String no_norm_flag = if defined(LRAA_no_norm) && LRAA_no_norm then "--no_norm" else ""
    String min_mapping_quality_flag = "--min_mapping_quality=" + select_first([LRAA_min_mapping_quality, 0])
    
    command <<<
        mkdir -p ~{OutDir}/Quant_noEM_minMapQ
    
        /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/Quant_noEM_minMapQ/LRAA.quant \
                                 --quant_only \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_full} \
                                 ~{min_mapping_quality_flag} --CPU 1
    >>>
    
    output {
        File lraaQuantExpr = "~{OutDir}/Quant_noEM_minMapQ/LRAA.quant.quant.expr"
        File lraaQuantTracking = "~{OutDir}/Quant_noEM_minMapQ/LRAA.quant.quant.tracking"
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
        Array[File] quantExprFiles
        Array[File] quantTrackingFiles
        String outputFilePrefix
        String docker
        Int memoryGB
        Int diskSizeGB
    }

    command <<<
        set -eo pipefail

        # Merge Quant Expression Files with header from the first file only
        quant_expr_output="~{outputFilePrefix}_merged_quant.expr"
        for file in ~{sep=' ' quantExprFiles}; do
            if [[ ! -f "$quant_expr_output" ]]; then
                # Copy the first file with header
                cp "$file" "$quant_expr_output"
            else
                # Append the rest without the header
                tail -n +2 "$file" >> "$quant_expr_output"
            fi
        done

        # Merge Quant Tracking Files with header from the first file only
        quant_tracking_output="~{outputFilePrefix}_merged_quant.tracking"
        for file in ~{sep=' ' quantTrackingFiles}; do
            if [[ ! -f "$quant_tracking_output" ]]; then
                # Copy the first file with header
                cp "$file" "$quant_tracking_output"
            else
                # Append the rest without the header
                tail -n +2 "$file" >> "$quant_tracking_output"
            fi
        done
    >>>

    output {
        File mergedQuantExprFile = "~{outputFilePrefix}_merged_quant.expr"
        File mergedQuantTrackingFile = "~{outputFilePrefix}_merged_quant.tracking"
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
        File inputBAM
        File referenceGenome
        Int? LRAA_min_mapping_quality
        Boolean? LRAA_no_norm
        Int cpu = 2
        Int numThreads = 4
        Int memoryGB = 32
        Int diskSizeGB = 1024
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        File referenceAnnotation_full
        String main_chromosomes = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
    }

    String OutDir = "LRAA_out"

    call splitBAMByChromosome {
        input:
            inputBAM = inputBAM,
            main_chromosomes = main_chromosomes,
            docker = docker,
            threads = numThreads,
            referenceGenome = referenceGenome,
            referenceAnnotation_full = referenceAnnotation_full,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }

    scatter (i in range(length(splitBAMByChromosome.chromosomeBAMs))) {
        call lraaPerChromosome {
            input:
                inputBAM = splitBAMByChromosome.chromosomeBAMs[i],
                referenceGenome = splitBAMByChromosome.chromosomeFASTAs[i],
                OutDir = OutDir,
                docker = docker,
                numThreads = numThreads,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality,
                referenceAnnotation_full = splitBAMByChromosome.fullAnnotations[i],
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
        }
    }
    call mergeResults {
        input:
            quantExprFiles = lraaPerChromosome.lraaQuantExpr,
            quantTrackingFiles = lraaPerChromosome.lraaQuantTracking,
            outputFilePrefix = "merged",
            docker = docker,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }
    
    output {
        File mergedQuantExpr = mergeResults.mergedQuantExprFile
        File mergedQuantTracking = mergeResults.mergedQuantTrackingFile
    }
}
