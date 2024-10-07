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
    String dockerImage = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    String outputDir = "LRAA_out"

    if (mode == "ID_ref_free_Quant_mode") {

        call IDRefFree.lraaWorkflow as IDRefFreeWorkflow {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = dockerImage,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm
        }

        call Quant.lraaWorkflow as QuantFreeWorkflow {
            input:
                inputBAMArray = IDRefFreeWorkflow.splitBAMs,
                referenceGenomeArray = IDRefFreeWorkflow.splitFASTAs,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = dockerImage,
                referenceAnnotation_full = IDRefFreeWorkflow.mergedReffreeGTF,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality
        }

        call Filtering.TranscriptFiltering as LRAA_ID_filtering_FreeWorkflow {
            input:
                gtf_path = IDRefFreeWorkflow.mergedReffreeGTF,
                expr_file_path = QuantFreeWorkflow.mergedQuantExpr,
                referenceGenome = referenceGenome,
                threshold = 1.0,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/filtertranscripts:latest"
        }

        call Quant.lraaWorkflow as QuantFree2Workflow {
            input:
                inputBAMArray = IDRefFreeWorkflow.splitBAMs,
                referenceGenomeArray = IDRefFreeWorkflow.splitFASTAs,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = dockerImage,
                referenceAnnotation_full = LRAA_ID_filtering_FreeWorkflow.filtered_gtf,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality
        }
    }

    if (mode == "ID_ref_guided_Quant_mode") {

        File guaranteedRefGuided = select_first([referenceGTF])

        call IDRefGuided.lraaWorkflow as IDRefGuidedWorkflow {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                referenceAnnotation_reduced = guaranteedRefGuided,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = dockerImage,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm
        }

        call Quant.lraaWorkflow as QuantGuidedWorkflow {
            input:
                inputBAMArray = IDRefGuidedWorkflow.splitBAMs,
                referenceGenomeArray = IDRefGuidedWorkflow.splitFASTAs,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = dockerImage,
                referenceAnnotation_full = IDRefGuidedWorkflow.mergedReffreeGTF,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality
        }

        call Filtering.TranscriptFiltering as LRAA_ID_filtering_GuidedWorkflow {
            input:
                gtf_path = IDRefGuidedWorkflow.mergedReffreeGTF,
                expr_file_path = QuantGuidedWorkflow.mergedQuantExpr,
                referenceGenome = referenceGenome,
                threshold = 1.0,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/filtertranscripts:latest"
        }

        call Quant.lraaWorkflow as QuantGuided2Workflow {
            input:
                inputBAMArray = IDRefGuidedWorkflow.splitBAMs,
                referenceGenomeArray = IDRefGuidedWorkflow.splitFASTAs,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = dockerImage,
                referenceAnnotation_full = LRAA_ID_filtering_GuidedWorkflow.filtered_gtf,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality
        }
    }

    if (mode == "Quant_only_mode") {

        File guaranteedRefQuantOnly = select_first([referenceGTF])

        call Quant.lraaWorkflow as QuantOnlyWorkflow {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = dockerImage,
                referenceAnnotation_full = guaranteedRefQuantOnly,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality
        }
    }

    output {
        File? UpdatedGTF = if (mode == "ID_ref_free_Quant_mode") then LRAA_ID_filtering_FreeWorkflow.filtered_gtf else if (mode == "ID_ref_guided_Quant_mode") then LRAA_ID_filtering_GuidedWorkflow.filtered_gtf else "null"
        File? Quant = if (mode == "ID_ref_free_Quant_mode") then QuantFree2Workflow.mergedQuantExpr else if (mode == "ID_ref_guided_Quant_mode") then QuantGuided2Workflow.mergedQuantExpr else QuantOnlyWorkflow.mergedQuantExpr
        File? Tracking = if (mode == "ID_ref_free_Quant_mode") then QuantFree2Workflow.mergedQuantTracking else if (mode == "ID_ref_guided_Quant_mode") then QuantGuided2Workflow.mergedQuantTracking else QuantOnlyWorkflow.mergedQuantTracking
    }
}
