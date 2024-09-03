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
        String IDQuantQuant_or_QuantOnly_or_Both
        Boolean? LRAA_no_norm
        Int? LRAA_min_mapping_quality
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        Int memoryGB
        Int diskSizeGB
    }

    String chrName = basename(inputBAM, '.bam')
    String no_norm_flag = if defined(LRAA_no_norm) && LRAA_no_norm then "--no_norm" else ""
    String min_mapping_quality_flag = "--min_mapping_quality=" + select_first([LRAA_min_mapping_quality, 0])
    
    command <<<
        mkdir -p ~{OutDir}/ID_reffree
        mkdir -p ~{OutDir}/ID_reduced
        mkdir -p ~{OutDir}/Quant_noEM_minMapQ
    
        # Use contig_names in the LRAA command
        if [[ ("~{IDQuantQuant_or_QuantOnly_or_Both}" == "IDQuantQuant" || "~{IDQuantQuant_or_QuantOnly_or_Both}" == "Both") && -z "~{referenceAnnotation_reduced}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/ID_reffree/LRAA_reffree \
                                     ~{no_norm_flag} --CPU 1
        fi
    
        if [[ ("~{IDQuantQuant_or_QuantOnly_or_Both}" == "IDQuantQuant" || "~{IDQuantQuant_or_QuantOnly_or_Both}" == "Both") && -f "~{referenceAnnotation_reduced}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/ID_reduced/LRAA_reduced \
                                     ~{no_norm_flag} \
                                     --gtf ~{referenceAnnotation_reduced} --CPU 1
        fi
    
        if [[ ("~{IDQuantQuant_or_QuantOnly_or_Both}" == "QuantOnly" || "~{IDQuantQuant_or_QuantOnly_or_Both}" == "Both") && -f "~{referenceAnnotation_full}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/Quant_noEM_minMapQ/LRAA.quant \
                                     --quant_only \
                                     ~{no_norm_flag} \
                                     --gtf ~{referenceAnnotation_full} \
                                     ~{min_mapping_quality_flag} --CPU 1
        fi
    >>>
    
    output {
        File lraaID_reffree_GTF = select_first(glob("~{OutDir}/ID_reffree/*_reffree.gtf"))
        File lraaID_reduced_GTF = select_first(glob("~{OutDir}/ID_reduced/*_reduced.gtf"))
        File lraaQuantExpr = select_first(glob("~{OutDir}/Quant_noEM_minMapQ/*.expr"))
        File lraaQuantTracking = select_first(glob("~{OutDir}/Quant_noEM_minMapQ/*.tracking"))
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
        Array[File] inputFiles
        String outputFile
        String docker
        Boolean isGTF
        Boolean isTracking
        Int memoryGB
        Int diskSizeGB
    }

    command <<<
        set -eo pipefail

        # Determine file extension based on file type
        ext=""
        if [[ ~{isGTF} == true ]]; then
            ext=".gtf"
        elif [[ ~{isTracking} == true ]]; then
            ext=".tracking"
        else
            ext=".expr"
        fi
        output_file="~{outputFile}"$ext
        touch $output_file

        # Write input files to a temporary file for better handling
        for file in ~{sep=" " inputFiles}; do
            echo $file >> input_files_list.txt
        done

        # Check if inputFiles array is not empty
        if [ ! -s input_files_list.txt ]; then
            echo "No input files provided."
            exit 1
        fi

        # Initialize a variable to track if the header has been added (for tracking or quant files)
        header_added=false

        # Merge files with conditions
        while IFS= read -r file; do
            if [[ -f "$file" ]]; then
                echo "Processing file: $file"
                if [[ ~{isGTF} == true ]]; then
                    # For GTF files, concatenate all files
                    cat $file >> $output_file
                else
                    # For tracking or quant files, handle header
                    if [[ $header_added == false ]]; then
                        cat $file >> $output_file
                        header_added=true
                    else
                        # Skip the header line(s) and append the rest
                        tail -n +2 $file >> $output_file
                    fi
                fi
            else
                echo "File $file does not exist."
            fi
        done < input_files_list.txt
    >>>

    output {
        File mergedFile = "~{outputFile}" + (if isGTF then ".gtf" else if isTracking then ".tracking" else ".expr")
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
        String IDQuantQuant_or_QuantOnly_or_Both
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
                IDQuantQuant_or_QuantOnly_or_Both = IDQuantQuant_or_QuantOnly_or_Both,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality,
                referenceAnnotation_reduced = select_first([splitBAMByChromosome.reducedAnnotations[i]]),
                referenceAnnotation_full = select_first([splitBAMByChromosome.fullAnnotations[i]]),
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
        }
    }


    call mergeResults as mergeReffreeGTF {
        input:
            inputFiles = lraaPerChromosome.lraaID_reffree_GTF,
            outputFile = "merged_reffree_ID",
            docker = docker,
            isGTF = true,
            isTracking = false,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }
    
    call mergeResults as mergeReducedGTF {
        input:
            inputFiles = lraaPerChromosome.lraaID_reduced_GTF,
            outputFile = "merged_reduced_ID",
            docker = docker,
            isGTF = true,
            isTracking = false,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }
    
    call mergeResults as mergeQuantExpr {
        input:
            inputFiles = lraaPerChromosome.lraaQuantExpr,
            outputFile = "merged_Quant",
            docker = docker,
            isGTF = false,
            isTracking = true,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }
    
    call mergeResults as mergeQuantTracking {
        input:
            inputFiles = lraaPerChromosome.lraaQuantTracking,
            outputFile = "merged_Quant.tracking",
            docker = docker,
            isGTF = false,
            isTracking = false,            
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
