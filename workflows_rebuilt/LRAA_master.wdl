version 1.0

import "LRAA_ID_ref_free.wdl" as IDRefFree
import "LRAA_ID_ref_guided.wdl" as IDRefGuided
import "LRAA_Quant.wdl" as Quant
import "LRAA_ID_filtering.wdl" as Filtering

workflow CombinedWorkflow {
    input {
        File inputBAM
        File referenceGenome
        File? referenceGTF
        String mode
        Int numThreads = 4
        Int memoryGB = 32
        String main_chromosomes = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
        Boolean? LRAA_no_norm
        Int? LRAA_min_mapping_quality = 0
    }

    Int diskSizeGB = 1024
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    String OutDir = "LRAA_out"

    if (mode == "ID_ref_free_Quant_mode") {

        call IDRefFree.lraaWorkflow as IDRefFree {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = docker,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm
        }

        call Quant.lraaWorkflow as QuantFree {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = docker,
                referenceAnnotation_full = IDRefFree.mergedReffreeGTF,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality
        }

        call Filtering.TranscriptFiltering as LRAA_ID_filtering_Free {
            input:
                gtf_path = IDRefFree.mergedReffreeGTF,
                expr_file_path = QuantFree.mergedQuantExpr,
                referenceGenome = referenceGenome,
                threshold = 1.0,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/filtertranscripts:latest"
        }

        call Quant.lraaWorkflow as QuantFree2 {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = docker,
                referenceAnnotation_full = LRAA_ID_filtering_Free.filtered_gtf,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality
        }
    }

    if (mode == "ID_ref_guided_Quant_mode") {

        if (defined(referenceGTF)) {
            call IDRefGuided.lraaWorkflow as IDRefGuided {
                input:
                    inputBAM = inputBAM,
                    referenceGenome = referenceGenome,
                    referenceAnnotation_reduced = referenceGTF,
                    numThreads = numThreads,
                    memoryGB = memoryGB,
                    diskSizeGB = diskSizeGB,
                    docker = docker,
                    main_chromosomes = main_chromosomes,
                    LRAA_no_norm = LRAA_no_norm
            }

            call Quant.lraaWorkflow as QuantGuided {
                input:
                    inputBAM = inputBAM,
                    referenceGenome = referenceGenome,
                    numThreads = numThreads,
                    memoryGB = memoryGB,
                    diskSizeGB = diskSizeGB,
                    docker = docker,
                    referenceAnnotation_full = IDRefGuided.mergedReducedGTF,
                    main_chromosomes = main_chromosomes,
                    LRAA_no_norm = LRAA_no_norm,
                    LRAA_min_mapping_quality = LRAA_min_mapping_quality
            }

            call Filtering.TranscriptFiltering as LRAA_ID_filtering_Guided {
                input:
                    gtf_path = IDRefGuided.mergedReducedGTF,
                    expr_file_path = QuantGuided.mergedQuantExpr,
                    referenceGenome = referenceGenome,
                    threshold = 1.0,
                    memoryGB = memoryGB,
                    diskSizeGB = diskSizeGB,
                    docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/filtertranscripts:latest"
            }

            call Quant.lraaWorkflow as QuantGuided2 {
                input:
                    inputBAM = inputBAM,
                    referenceGenome = referenceGenome,
                    numThreads = numThreads,
                    memoryGB = memoryGB,
                    diskSizeGB = diskSizeGB,
                    docker = docker,
                    referenceAnnotation_full = LRAA_ID_filtering_Guided.filtered_gtf,
                    main_chromosomes = main_chromosomes,
                    LRAA_no_norm = LRAA_no_norm,
                    LRAA_min_mapping_quality = LRAA_min_mapping_quality
            }
        } else {
            fail("referenceGTF must be provided for ID_ref_guided_Quant_mode")
        }
    }

    if (mode == "Quant_only") {

        if (defined(referenceGTF)) {
            call Quant.lraaWorkflow as QuantOnly {
                input:
                    inputBAM = inputBAM,
                    referenceGenome = referenceGenome,
                    numThreads = numThreads,
                    memoryGB = memoryGB,
                    diskSizeGB = diskSizeGB,
                    docker = docker,
                    referenceAnnotation_full = referenceGTF,
                    main_chromosomes = main_chromosomes,
                    LRAA_no_norm = LRAA_no_norm,
                    LRAA_min_mapping_quality = LRAA_min_mapping_quality
            }
        } else {
            fail("referenceGTF must be provided for Quant_only mode")
        }
    }

    output {
        File? UpdatedGTF = if (mode == "ID_ref_free_Quant_mode") then LRAA_ID_filtering_Free.filtered_gtf else if (mode == "ID_ref_guided_Quant_mode") then LRAA_ID_filtering_Guided.filtered_gtf else "null"
        File? Quant = if (mode == "ID_ref_free_Quant_mode") then QuantFree2.mergedQuantExpr else if (mode == "ID_ref_guided_Quant_mode") then QuantGuided2.mergedQuantExpr else QuantOnly.mergedQuantExpr
        File? Tracking = if (mode == "ID_ref_free_Quant_mode") then QuantFree2.mergedQuantTracking else if (mode == "ID_ref_guided_Quant_mode") then QuantGuided2.mergedQuantTracking else QuantOnly.mergedQuantTracking
    }
}
