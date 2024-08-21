version 1.0

task splitBAMByChromosome {
    input {
        File inputBAM
        String docker = "biocontainers/samtools:v1.9-4-deb_cv1"
        Array[String] chromosomes
    }
    command <<<
        mkdir -p split_bams
        for chr in ~{sep=' ' chromosomes}; do
            samtools view -b ~{inputBAM} $chr > split_bams/$chr.bam
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
        String? LRAA_min_mapping_quality_flag
        String? no_norm_flag
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        File monitoringScript
    }
    command <<<
        bash ~{monitoringScript} > monitoring.log &
        
        mkdir -p ~{OutDir}/ID ~{OutDir}/ID_reduced ~{OutDir}/Quant ~{OutDir}/Quant_noEM ~{OutDir}/Quant_minMapQ ~{OutDir}/Quant_noEM_minMapQ

        if [[ "~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/ID/LRAA \
                                     ~{no_norm_flag} --CPU ~{numThreads}
        fi

        if [[ ("~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both") && -n "~{referenceAnnotation_reduced}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/ID_reduced/LRAA_reduced \
                                     ~{no_norm_flag} \
                                     --gtf ~{referenceAnnotation_reduced} --CPU ~{numThreads}
        fi

        if [[ ("~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both") && -n "~{referenceAnnotation_full}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                     --bam ~{inputBAM} \
                                     --output_prefix ~{OutDir}/Quant_noEM_minMapQ/LRAA.noEM.minMapQ \
                                     --quant_only \
                                     ~{no_norm_flag} \
                                     --gtf ~{referenceAnnotation_full} \
                                     ~{LRAA_min_mapping_quality_flag} --CPU ~{numThreads}
        fi
    >>>
    output {
        File? lraaIDGTF = glob("~{OutDir}/ID/*.gtf")[0]
        File? lraaIDReducedGTF = glob("~{OutDir}/ID_reduced/*.gtf")[0]
        File? lraaQuantExpr = glob("~{OutDir}/Quant_noEM_minMapQ/*.quant.expr")[0]
        File? lraaQuantTracking = glob("~{OutDir}/Quant_noEM_minMapQ/*.quant.tracking")[0]
    }
    runtime {
        docker: docker
    }
}

task mergeResults {
    input {
        Array[File] inputFiles
        String outputFile
        String docker = "ubuntu:18.04"
    }
    command <<<
        if [ ! -z "~{inputFiles}" ]; then
            # Take the header from the first file
            head -n 1 ~{inputFiles[0]} > ~{outputFile}
            for file in ~{sep=' ' inputFiles}; do
                # Skip the header for subsequent files
                tail -n +2 $file >> ~{outputFile}
            done
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
        Array[String] chromosomes
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
        File monitoringScript
    }

    String OutDir = "LRAA_out"
    String LRAA_min_mapping_quality_flag = if (defined(LRAA_min_mapping_quality)) then "--min_mapping_quality=" + LRAA_min_mapping_quality else ""
    String no_norm_flag = if (defined(LRAA_no_norm) && LRAA_no_norm) then "--no_norm" else ""

    call splitBAMByChromosome {
        input:
            inputBAM = inputBAM,
            chromosomes = chromosomes,
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
    Array[File] idGTFFiles = flatten(select_all(lraaPerChromosome.lraaIDGTF))
    Array[File] idReducedGTFFiles = flatten(select_all(lraaPerChromosome.lraaIDReducedGTF))

    call mergeResults as mergeIDGTF {
        input:
            inputFiles = idGTFFiles,
            outputFile = OutDir + "/merged_ID.gtf",
            docker = docker
    }

    call mergeResults as mergeIDReducedGTF {
        input:
            inputFiles = idReducedGTFFiles,
            outputFile = OutDir + "/merged_ID_reduced.gtf",
            docker = docker
    }

    # Merge Quant results (.expr and .tracking files)
    Array[File] quantExprFiles = select_all(lraaPerChromosome.lraaQuantExpr)
    Array[File] quantTrackingFiles = select_all(lraaPerChromosome.lraaQuantTracking)

    call mergeResults as mergeQuantExpr {
        input:
            inputFiles = quantExprFiles,
            outputFile = OutDir + "/merged_Quant.expr",
            docker = docker
    }

    call mergeResults as mergeQuantTracking {
        input:
            inputFiles = quantTrackingFiles,
            outputFile = OutDir + "/merged_Quant.tracking",
            docker = docker
    }

    output {
        File mergedIDGTF = mergeIDGTF.mergedFile
        File mergedIDReducedGTF = mergeIDReducedGTF.mergedFile
        File mergedQuantExpr = mergeQuantExpr.mergedFile
        File mergedQuantTracking = mergeQuantTracking.mergedFile
    }
}