version 1.0

task splitBAMByChromosome {
    input {
        File inputBAM
        String main_chromosomes
        String docker
        Int threads
    }

    command <<<
        set -eo pipefail
        mkdir split_bams
        samtools view -@ ~{threads} -H ~{inputBAM} > header.sam
        for chr in ~{main_chromosomes}; do
            samtools view -@ ~{threads} -b ~{inputBAM} $chr > split_bams/$chr.bam
        done
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
    }

    command <<<
        set -eo pipefail
        mkdir -p ~{OutDir}
        lraa \
            --bam ~{inputBAM} \
            --genome ~{referenceGenome} \
            --out ~{OutDir} \
            --threads ~{numThreads} \
            --mode ~{ID_or_Quant_or_Both} \
            ~{"--no-norm" if defined(LRAA_no_norm) and LRAA_no_norm == true else ""} \
            ~{"--min-mapping-quality " + LRAA_min_mapping_quality if defined(LRAA_min_mapping_quality) else ""} \
            ~{"--annotation-reduced " + referenceAnnotation_reduced if defined(referenceAnnotation_reduced) else ""} \
            ~{"--annotation-full " + referenceAnnotation_full if defined(referenceAnnotation_full) else ""}
    >>>

    output {
        File? lraaID_reffree_GTF = glob("~{OutDir}/*_reffree.gtf")[0]
        File? lraaID_reduced_GTF = glob("~{OutDir}/*_reduced.gtf")[0]
        File? lraaQuantExpr = glob("~{OutDir}/*.expr")[0]
        File? lraaQuantTracking = glob("~{OutDir}/*.tracking")[0]
    }

    runtime {
        docker: docker
    }
}

task mergeResults {
    input {
        Array[File?] inputFiles
        String outputFile
        String docker
        Boolean isGTF
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
        File mergedFile = "~{outputFile}" + ext
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
            outputFile = OutDir + "/merged_reffree_ID",
            docker = docker,
            isGTF = true
    }

    # Collect and merge reduced GTF files
    Array[File?] reducedGTFFiles = select_all(lraaPerChromosome.lraaID_reduced_GTF)
    call mergeResults as mergeReducedGTF {
        input:
            inputFiles = reducedGTFFiles,
            outputFile = OutDir + "/merged_reduced_ID",
            docker = docker,
            isGTF = true
    }

    Array[File?] quantExprFiles = select_all(lraaPerChromosome.lraaQuantExpr)
    Array[File?] quantTrackingFiles = select_all(lraaPerChromosome.lraaQuantTracking)

    call mergeResults as mergeQuantExpr {
        input:
            inputFiles = quantExprFiles,
            outputFile = OutDir + "/merged_Quant",
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
