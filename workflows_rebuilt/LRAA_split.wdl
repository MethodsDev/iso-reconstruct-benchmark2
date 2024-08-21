version 1.0

task splitBAMByChromosome {
    input {
        File inputBAM
        String main_chromosomes = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        Int threads
    }
    command <<<
        set -eo pipefail

        if [ ! -f "~{inputBAM}.bai" ]; then
            samtools index -@ ~{threads} ~{inputBAM}
        fi

        chromosomes=$(samtools idxstats ~{inputBAM} | cut -f1 | grep -vE '^$|chrM')

        mkdir -p split_bams
        other_contigs_bam="split_bams/other_contigs.bam"
        touch $other_contigs_bam

        for chr in $chromosomes; do
            if [[ " ~{main_chromosomes} " =~ .*\ $chr\ .* ]]; then
                samtools view -@ ~{threads} -b ~{inputBAM} $chr > split_bams/$chr.bam
            else
                samtools view -@ ~{threads} -b ~{inputBAM} $chr >> $other_contigs_bam
            fi
        done

        samtools index -@ ~{threads} $other_contigs_bam
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
        Boolean? LRAA_no_norm
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        Int? LRAA_min_mapping_quality
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    }
    command <<<
        ${if defined(LRAA_no_norm) && LRAA_no_norm then "--no_norm" else ""}
        ${if defined(LRAA_min_mapping_quality) then "--min_mapping_quality=" + LRAA_min_mapping_quality else ""}

        mkdir -p ~{OutDir}/ID_reffree
        mkdir -p ~{OutDir}/ID_reduced
        mkdir -p ~{OutDir}/Quant_noEM_minMapQ

        if [[ "~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/ID_reffree/LRAA \
                                 ${if defined(LRAA_no_norm) && LRAA_no_norm then "--no_norm" else ""} --CPU ~{numThreads}
        fi

        if [[ ("~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both") && -n "~{referenceAnnotation_reduced}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/ID_reduced/LRAA_reduced \
                                 ${if defined(LRAA_no_norm) && LRAA_no_norm then "--no_norm" else ""} \
                                 --gtf ~{referenceAnnotation_reduced} --CPU ~{numThreads} 
        fi

        if [[ ("~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both") && -n "~{referenceAnnotation_full}" && -n "~{LRAA_min_mapping_quality}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/Quant_noEM_minMapQ/LRAA.noEM.minMapQ \
                                 --quant_only \
                                 ${if defined(LRAA_no_norm) && LRAA_no_norm then "--no_norm" else ""} \
                                 --gtf ~{referenceAnnotation_full} \
                                 ${if defined(LRAA_min_mapping_quality) then "--min_mapping_quality=" + LRAA_min_mapping_quality else ""} --CPU ~{numThreads}
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
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 2048
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
                referenceAnnotation_reduced = referenceAnnotation_reduced,
                referenceAnnotation_full = referenceAnnotation_full
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
