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

        call Quant.lraaWorkflow as Quant {
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

        call Filtering.TranscriptFiltering as LRAA_ID_filtering {
            input:
                gtf_path = IDRefFree.mergedReffreeGTF,
                expr_file_path = Quant.mergedQuantExpr,
                referenceGenome = referenceGenome,
                threshold = 1.0,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/filtertranscripts:latest"
        }

        call Quant.lraaWorkflow as Quant2 {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = docker,
                referenceAnnotation_full = LRAA_ID_filtering.filtered_gtf,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality
        }
    }

    if (mode == "ID_ref_guided_Quant_mode") {

        call IDRefGuided.lraaWorkflow as IDRefGuided {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                referenceGTF = referenceGTF,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = docker,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm
        }

        call Quant.lraaWorkflow as Quant {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = docker,
                referenceAnnotation_full = IDRefGuided.mergedRefguidedGTF,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality
        }

        call Filtering.TranscriptFiltering as LRAA_ID_filtering {
            input:
                gtf_path = IDRefGuided.mergedRefguidedGTF,
                expr_file_path = Quant.mergedQuantExpr,
                referenceGenome = referenceGenome,
                threshold = 1.0,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/filtertranscripts:latest"
        }

        call Quant.lraaWorkflow as Quant2 {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = docker,
                referenceAnnotation_full = LRAA_ID_filtering.filtered_gtf,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality
        }
    }

    output {
        File? UpdatedGTF = LRAA_ID_filtering.filtered_gtf
        File? Quant = Quant2.mergedQuantExpr
        File? Tracking = Quant2.mergedQuantTracking
    }
}
