version 1.0

task splitBAMByChromosome {
    input {
        File inputBAM
        String main_chromosomes = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
        String docker = "biocontainers/samtools:v1.9-4-deb_cv1"
    }
    command <<<
        set -eo pipefail

        # Index the BAM file if not already indexed
        if [ ! -f "~{inputBAM}.bai" ]; then
            samtools index ~{inputBAM}
        fi

        # Extract chromosome names
        chromosomes=$(samtools idxstats ~{inputBAM} | cut -f1 | grep -vE '^$|chrM')

        mkdir -p split_bams
        other_contigs_bam="split_bams/other_contigs.bam"
        touch $other_contigs_bam

        for chr in $chromosomes; do
            if [[ " ~{main_chromosomes} " =~ .*\ $chr\ .* ]]; then
                # Split main chromosomes
                samtools view -b ~{inputBAM} $chr > split_bams/$chr.bam
            else
                # Combine other contigs/scaffolds
                samtools view -b ~{inputBAM} $chr >> $other_contigs_bam
            fi
        done

        # Index the combined other contigs BAM file
        samtools index $other_contigs_bam
    >>>
    output {
        Array[File] chromosomeBAMs = glob("split_bams/*.bam")
    }
    runtime {
        docker: docker
    }
}

task lraaPerChromosome {
    input {
        String ID_or_Quant_or_Both
        File referenceGenome
        File inputBAM
        String OutDir
        Int numThreads
        Boolean no_norm_flag
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String? LRAA_min_mapping_quality_flag
    }

    # Define command line flags based on conditions before the command section
    String no_norm_cmd = if no_norm_flag then "--no_norm" else ""
    String min_mapping_quality_cmd = if defined(LRAA_min_mapping_quality_flag) then "--min_mapping_quality " + LRAA_min_mapping_quality_flag else ""

    command <<<
        mkdir -p ~{OutDir}/ID_reduced
        mkdir -p ~{OutDir}/Quant_noEM_minMapQ

        if [[ "~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            if [ -n "~{referenceAnnotation_reduced}" ]; then
                /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                         --bam ~{inputBAM} \
                                         --output_prefix ~{OutDir}/ID_reduced/LRAA_reduced \
                                         ~{no_norm_cmd} \
                                         --gtf ~{referenceAnnotation_reduced} --CPU ~{numThreads}
            fi
        fi

        if [[ "~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            if [ -n "~{referenceAnnotation_full}" ]; then
                /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                         --bam ~{inputBAM} \
                                         --output_prefix ~{OutDir}/Quant_noEM_minMapQ/LRAA.noEM.minMapQ \
                                         --quant_only \
                                         ~{no_norm_cmd} \
                                         --gtf ~{referenceAnnotation_full} \
                                         ~{min_mapping_quality_cmd} --CPU ~{numThreads}
            fi
        fi
    >>>

    output {
        File? lraaIDGTF = if (length(glob("~{OutDir}/ID/*.gtf")) > 0) then glob("~{OutDir}/ID/*.gtf")[0] else None
        File? lraaIDReducedGTF = if (length(glob("~{OutDir}/ID_reduced/*.gtf")) > 0) then glob("~{OutDir}/ID_reduced/*.gtf")[0] else None
        File? lraaQuantExpr = if (length(glob("~{OutDir}/Quant_noEM_minMapQ/*.quant.expr")) > 0) then glob("~{OutDir}/Quant_noEM_minMapQ/*.quant.expr")[0] else None
        File? lraaQuantTracking = if (length(glob("~{OutDir}/Quant_noEM_minMapQ/*.quant.tracking")) > 0) then glob("~{OutDir}/Quant_noEM_minMapQ/*.quant.tracking")[0] else None
    }
    runtime {
        docker: docker
    }
}

task mergeResults {
    input {
        Array[File?] inputFiles
        String outputFile
        String docker = "ubuntu:18.04"
        Boolean isGTF = false
    }
    command <<<
        if [ ! -z "~{inputFiles}" ]; then
            # Filter out nulls and take the header from the first file
            filtered_files=(~{sep=' ' inputFiles})
            if [[ "~{isGTF}" == "true" ]]; then
                cat "${filtered_files[@]}" > ~{outputFile}
            else
                head -n 1 "${filtered_files[0]}" > ~{outputFile}
                for file in "${filtered_files[@]}"; do
                    # Skip the header for subsequent files
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
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 2048
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
        String main_chromosomes = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
    }

    String OutDir = "LRAA_out"
    String LRAA_min_mapping_quality_flag = if (defined(LRAA_min_mapping_quality)) then "--min_mapping_quality=" + LRAA_min_mapping_quality else ""
    String no_norm_flag = if (defined(LRAA_no_norm) && LRAA_no_norm) then "--no_norm" else ""

    call splitBAMByChromosome {
        input:
            inputBAM = inputBAM,
            main_chromosomes = main_chromosomes,
            docker = docker
    }

    scatter (chrBAM in splitBAMByChromosome.chromosomeBAMs) {
        call lraaPerChromosome {
            input:
                inputBAM = chrBAM,
                referenceGenome = referenceGenome,
                OutDir = OutDir,
                docker = docker,
                numThreads = numThreads,
                ID_or_Quant_or_Both = ID_or_Quant_or_Both,
                LRAA_min_mapping_quality_flag = LRAA_min_mapping_quality_flag,
                no_norm_flag = no_norm_flag,
                referenceAnnotation_reduced = referenceAnnotation_reduced,
                referenceAnnotation_full = referenceAnnotation_full,
                monitoringScript = monitoringScript
        }
    }

    # Merge ID results (GTF files)
    Array[File?] idGTFFiles = select_all(lraaPerChromosome.lraaIDGTF)
    Array[File?] idReducedGTFFiles = select_all(lraaPerChromosome.lraaIDReducedGTF)
    
    call mergeResults as mergeIDGTF {
        input:
            inputFiles = idGTFFiles,
            outputFile = OutDir + "/merged_ID.gtf",
            docker = docker,
            isGTF = true
    }
    
    call mergeResults as mergeIDReducedGTF {
        input:
            inputFiles = idReducedGTFFiles,
            outputFile = OutDir + "/merged_ID_reduced.gtf",
            docker = docker,
            isGTF = true
    }

    # Merge Quant results (.expr and .tracking files)
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
        File mergedIDGTF = mergeIDGTF.mergedFile
        File mergedIDReducedGTF = mergeIDReducedGTF.mergedFile
        File mergedQuantExpr = mergeQuantExpr.mergedFile
        File mergedQuantTracking = mergeQuantTracking.mergedFile
    }
}
