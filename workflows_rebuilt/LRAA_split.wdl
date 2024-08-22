version 1.0

task splitBAMByChromosome {
    input {
        File inputBAM
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String main_chromosomes = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        Int threads
    }
    command <<<
        set -eo pipefail
        
        if [ ! -f "~{inputBAM}.bai" ]; then
            samtools index -@ ~{threads} ~{inputBAM}
        fi
        
        mkdir -p split_bams
        mkdir -p split_gtf_reduced
        mkdir -p split_gtf_full
        
        main_chromosomes="~{main_chromosomes}"
        
        for chr in $main_chromosomes; do
            samtools view -@ ~{threads} -b ~{inputBAM} $chr > split_bams/$chr.bam
            if [ -f "~{referenceAnnotation_reduced}" ]; then
                grep "^$chr" ~{referenceAnnotation_reduced} > split_gtf_reduced/$chr.gtf || echo "" > split_gtf_reduced/$chr.gtf
            fi
            if [ -f "~{referenceAnnotation_full}" ]; then
                grep "^$chr" ~{referenceAnnotation_full} > split_gtf_full/$chr.gtf || echo "" > split_gtf_full/$chr.gtf
            fi
        done
        
        # Handle other contigs by first creating a list of all contigs in the BAM file
        samtools idxstats ${inputBAM} | cut -f1 | grep -v -E $(echo $main_chromosomes | sed 's/ /|/g') > other_contigs.txt
        
        # Then, use samtools view to create a BAM file for other contigs
        samtools view -@ ${threads} -b ${inputBAM} -o split_bams/other_contigs.bam -L other_contigs.txt

        if [ -f "~{referenceAnnotation_reduced}" ]; then
            grep -vE "($(echo $main_chromosomes | sed 's/ /|/g'))" ~{referenceAnnotation_reduced} > split_gtf_reduced/other_contigs.gtf
        fi
        if [ -f "~{referenceAnnotation_full}" ]; then
            grep -vE "($(echo $main_chromosomes | sed 's/ /|/g'))" ~{referenceAnnotation_full} > split_gtf_full/other_contigs.gtf
        fi
    >>>
    output {
        Array[File] chromosomeBAMs = glob("split_bams/*.bam")
        Array[File] chromosomeGTFs_reduced = glob("split_gtf_reduced/*.gtf")
        Array[File] chromosomeGTFs_full = glob("split_gtf_full/*.gtf")
    }
    runtime {
        docker: docker
    }
}

task FilterGTF {
    input {
        Array[File] gtfFiles
        String chrName
    }

    command <<<
        #!/bin/bash
        for gtf in "~{sep=' ' gtfFiles}"; do
            if [[ $(basename "$gtf" .gtf) == "~{chrName}" ]]; then
                echo "$gtf"
                break
            fi
        done
    >>>

    output {
        File? selectedGTF = read_string(stdout())
    }

    runtime {
        docker: "ubuntu:latest"
    }
}


task lraaPerChromosome {
    input {
        String ID_or_Quant_or_Both
        File referenceGenome
        File inputBAM
        String OutDir
        Int numThreads
        Boolean? LRAA_no_norm
        Array[File] referenceAnnotation_reduced_chroms
        Array[File] referenceAnnotation_full_chroms
        Int? LRAA_min_mapping_quality
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        File? FilterReducedGTF
        File? FilterFullGTF
    }

    String chrName = basename(inputBAM, '.bam')
    String no_norm_flag = if defined(LRAA_no_norm) && LRAA_no_norm then "--no_norm" else ""
    String min_mapping_quality_flag = if defined(LRAA_min_mapping_quality) then "--min_mapping_quality=" + LRAA_min_mapping_quality else ""

    command <<<
        mkdir -p ~{OutDir}/ID_reffree
        mkdir -p ~{OutDir}/ID_reduced
        mkdir -p ~{OutDir}/Quant_noEM_minMapQ
    
        # Extract chromosome name or handle 'other_contigs'
        contig_names=$(if [[ "~{inputBAM}" =~ .*other_contigs.bam ]]; then
            awk '{print $1}' ~{FilterFullGTF} | sort | uniq | tr '\n' ','
        else
            echo "~{chrName}"
        fi)
        contig_names=${contig_names%,} # Remove trailing comma
    
        # Use contig_names in the LRAA command
        if [[ ("~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both") && -f "~{FilterReducedGTF}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/ID_reffree/LRAA \
                                     ~{no_norm_flag} --CPU ~{numThreads} \
                                     --contig $contig_names
        fi
    
        if [[ ("~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both") && -f "~{FilterReducedGTF}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/ID_reduced/LRAA_reduced \
                                     ~{no_norm_flag} \
                                     --gtf ~{FilterReducedGTF} --CPU ~{numThreads} \
                                     --contig $contig_names
        fi
    
        if [[ ("~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both") && -f "~{FilterFullGTF}" && defined(LRAA_min_mapping_quality) ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/Quant_noEM_minMapQ/LRAA.noEM.minMapQ \
                                     --quant_only \
                                     ~{no_norm_flag} \
                                     --gtf ~{FilterFullGTF} \
                                     ~{min_mapping_quality_flag} --CPU ~{numThreads} \
                                     --contig $contig_names
        fi
    >>>

    output {
        File? lraaID_reffree_GTF = if (length(glob("~{OutDir}/ID_reffree/*.gtf")) > 0) then glob("~{OutDir}/ID_reffree/*.gtf")[0] else ""
        File? lraaID_reduced_GTF = if (length(glob("~{OutDir}/ID_reduced/*.gtf")) > 0) then glob("~{OutDir}/ID_reduced/*.gtf")[0] else ""
        File? lraaQuantExpr = if (length(glob("~{OutDir}/Quant_noEM_minMapQ/*.quant.expr")) > 0) then glob("~{OutDir}/Quant_noEM_minMapQ/*.quant.expr")[0] else ""
        File? lraaQuantTracking = if (length(glob("~{OutDir}/Quant_noEM_minMapQ/*.quant.tracking")) > 0) then glob("~{OutDir}/Quant_noEM_minMapQ/*.quant.tracking")[0] else ""
    }

    runtime {
        docker: docker
    }
}

task mergeResults {
    input {
        Array[File?] inputFiles
        String outputFile
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        Boolean isGTF = false
    }
    command <<<
        if [ ! -z "~{inputFiles}" ]; then
            filtered_files=(~{sep=' ' inputFiles})
            if [[ "~{isGTF}" == "true" ]]; then
                cat "${filtered_files[@]}" > ~{outputFile}
            else
                head -n 1 "${filtered_files[0]}" > ~{outputFile}
                for file in "${filtered_files[@]}"; do
                    tail -n +2 $file >> ~{outputFile}
                done
            fi
        fi
    >>>
    output {
        File mergedFile = outputFile
    }
    runtime {
        docker: docker
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
        Int diskSizeGB = 512
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
            threads = numThreads
    }

    scatter (chrBAM in splitBAMByChromosome.chromosomeBAMs) {
        String chrName = basename(chrBAM, '.bam')
    
        call FilterGTF as FilterReducedGTF {
            input:
                gtfFiles = splitBAMByChromosome.chromosomeGTFs_reduced,
                chrName = chrName
        }
    
        call FilterGTF as FilterFullGTF {
            input:
                gtfFiles = splitBAMByChromosome.chromosomeGTFs_full,
                chrName = chrName
        }
    
        call lraaPerChromosome {
            input:
                inputBAM = chrBAM,
                referenceGenome = referenceGenome,
                OutDir = OutDir,
                docker = docker,
                numThreads = numThreads,
                ID_or_Quant_or_Both = ID_or_Quant_or_Both,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality,
                referenceAnnotation_reduced_chroms = select_all([FilterReducedGTF.selectedGTF]),
                referenceAnnotation_full_chroms = select_all([FilterFullGTF.selectedGTF])
        }
    }


    # Collect and merge reffree GTF files
    Array[File?] reffreeGTFFiles = select_all(lraaPerChromosome.lraaID_reffree_GTF)
    call mergeResults as mergeReffreeGTF {
        input:
            inputFiles = reffreeGTFFiles,
            outputFile = OutDir + "/merged_reffree_ID.gtf",
            docker = docker,
            isGTF = true
    }

    # Collect and merge reduced GTF files
    Array[File?] reducedGTFFiles = select_all(lraaPerChromosome.lraaID_reduced_GTF)
    call mergeResults as mergeReducedGTF {
        input:
            inputFiles = reducedGTFFiles,
            outputFile = OutDir + "/merged_reduced_ID.gtf",
            docker = docker,
            isGTF = true
    }

    Array[File?] quantExprFiles = select_all(lraaPerChromosome.lraaQuantExpr)
    Array[File?] quantTrackingFiles = select_all(lraaPerChromosome.lraaQuantTracking)

    call mergeResults as mergeQuantExpr {
        input:
            inputFiles = quantExprFiles,
            outputFile = OutDir + "/merged_Quant.expr",
            docker = docker,
            isGTF = false
    }

    call mergeResults as mergeQuantTracking {
        input:
            inputFiles = quantTrackingFiles,
            outputFile = OutDir + "/merged_Quant.tracking",
            docker = docker,
            isGTF = false
    }


    output {
        File mergedReffreeGTF = mergeReffreeGTF.mergedFile
        File mergedReducedGTF = mergeReducedGTF.mergedFile
        File mergedQuantExpr = mergeQuantExpr.mergedFile
        File mergedQuantTracking = mergeQuantTracking.mergedFile
    }
}
