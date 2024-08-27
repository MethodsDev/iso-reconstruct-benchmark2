version 1.0

task splitBAMByChromosome {
    input {
        File inputBAM
        String main_chromosomes
        String docker
        Int threads
        File referenceGenome
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
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
            
        # Generate chromosome-specific GTF for reduced annotation, if available
        if [ -f "~{referenceAnnotation_reduced}" ]; then
            cat ~{referenceAnnotation_reduced} | perl -lane 'if ($F[0] eq "'$chr'") { print; }' > split_bams/$chr.reduced.annot.gtf
        fi
        
        # Generate chromosome-specific GTF for full annotation, if available
        if [ -f "~{referenceAnnotation_full}" ]; then
            cat ~{referenceAnnotation_full} | perl -lane 'if ($F[0] eq "'$chr'") { print; }' > split_bams/$chr.full.annot.gtf
        fi
        done
    >>>

    output {
        Array[File] chromosomeBAMs = glob("split_bams/*.bam")
        Array[File] chromosomeFASTAs = glob("split_bams/*.genome.fasta")
        Array[File] reducedAnnotations = glob("split_bams/*.reduced.annot.gtf")
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
        String ID_or_Quant_or_Both
        Boolean? LRAA_no_norm
        Int? LRAA_min_mapping_quality
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        Int memoryGB
        Int diskSizeGB
    }

    String chrName = basename(inputBAM, '.bam')
    String no_norm_flag = if defined(LRAA_no_norm) && LRAA_no_norm then "--no_norm" else ""
    String min_mapping_quality_flag = if defined(LRAA_min_mapping_quality) then "--min_mapping_quality=" + LRAA_min_mapping_quality else ""
    
    command <<<
        mkdir -p ~{OutDir}/ID_reffree
        mkdir -p ~{OutDir}/ID_reduced
        mkdir -p ~{OutDir}/Quant_noEM_minMapQ
    
        # Use contig_names in the LRAA command
        if [[ ("~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both") && -z "~{referenceAnnotation_reduced}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/ID_reffree/LRAA_reffree \
                                     ~{no_norm_flag} --CPU 1
        fi
    
        if [[ ("~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both") && -f "~{referenceAnnotation_reduced}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/ID_reduced/LRAA_reduced \
                                     ~{no_norm_flag} \
                                     --gtf ~{referenceAnnotation_reduced} --CPU 1
        fi
    
        if [[ ("~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both") && -f "~{referenceAnnotation_full}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/Quant_noEM_minMapQ/LRAA.noEM.minMapQ \
                                     --quant_only \
                                     ~{no_norm_flag} \
                                     --gtf ~{referenceAnnotation_full} \
                                     ~{min_mapping_quality_flag} --CPU 1
        fi
    >>>

    output {
        Array[File?] lraaID_reffree_GTF = select_first(glob("~{OutDir}/ID_reffree/*_reffree.gtf"))
        Array[File?] lraaID_reduced_GTF = select_first(glob("~{OutDir}/ID_reduced/*_reduced.gtf"))
        Array[File?] lraaQuantExpr = select_first(glob("~{OutDir}/Quant_noEM_minMapQ/*.expr"))
        Array[File?] lraaQuantTracking = select_first(glob("~{OutDir}/Quant_noEM_minMapQ/*.tracking"))
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
        Array[File?] inputFiles
        String outputFile
        String docker
        Boolean isGTF
        Int memoryGB
        Int diskSizeGB
    }

    command <<<
        set -eo pipefail

        # Determine file extension based on file type
        ext=""
        if [[ ~{isGTF} == true ]]; then
            ext=".gtf"
        else
            ext=".expr"
        fi

        # Merge files
        for file in ~{sep=" " inputFiles}; do
            if [[ -f "$file" ]]; then
                cat $file >> ~{outputFile}$ext
            fi
        done
    >>>

    output {
        File mergedFile = "~{outputFile}" + (if isGTF then ".gtf" else ".expr")
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
        String ID_or_Quant_or_Both
        Int? LRAA_min_mapping_quality
        Boolean? LRAA_no_norm
        Int cpu = 2
        Int numThreads = 4
        Int memoryGB = 32
        Int diskSizeGB = 1024
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
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
            referenceAnnotation_reduced = referenceAnnotation_reduced,
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
                ID_or_Quant_or_Both = ID_or_Quant_or_Both,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality,
                referenceAnnotation_reduced = select_first([splitBAMByChromosome.reducedAnnotations[i]]),
                referenceAnnotation_full = select_first([splitBAMByChromosome.fullAnnotations[i]]),
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
        }
    }

    # Collect and merge reffree GTF files
    Array[File?] reffreeGTFFiles = flatten(select_all(lraaPerChromosome.lraaID_reffree_GTF))
    call mergeResults as mergeReffreeGTF {
        input:
            inputFiles = reffreeGTFFiles,
            outputFile = OutDir + "/merged_reffree_ID",
            docker = docker,
            isGTF = true,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }

    # Collect and merge reduced GTF files
    Array[File?] reducedGTFFiles = flatten(select_all(lraaPerChromosome.lraaID_reduced_GTF))
    call mergeResults as mergeReducedGTF {
        input:
            inputFiles = reducedGTFFiles,
            outputFile = OutDir + "/merged_reduced_ID",
            docker = docker,
            isGTF = true,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }

    Array[File?] quantExprFiles = flatten(select_all(lraaPerChromosome.lraaQuantExpr))
    Array[File?] quantTrackingFiles = flatten(select_all(lraaPerChromosome.lraaQuantTracking))

    call mergeResults as mergeQuantExpr {
        input:
            inputFiles = quantExprFiles,
            outputFile = OutDir + "/merged_Quant",
            docker = docker,
            isGTF = false,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }

    call mergeResults as mergeQuantTracking {
        input:
            inputFiles = quantTrackingFiles,
            outputFile = OutDir + "/merged_Quant.tracking",
            docker = docker,
            isGTF = false,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }

    output {
        File mergedReffreeGTF = mergeReffreeGTF.mergedFile
        File mergedReducedGTF = mergeReducedGTF.mergedFile
        File mergedQuantExpr = mergeQuantExpr.mergedFile
        File mergedQuantTracking = mergeQuantTracking.mergedFile
    }
}
