version 1.0

import "ID_ref_free.wdl" as IDRefFree
import "ID_ref_guided.wdl" as IDRefGuided
import "Quant.wdl" as Quant
import "Filtering.wdl" as Filtering

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
        LRAA_min_mapping_quality? = 0
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


        call LRAA_ID_filtering.TranscriptFiltering as LRAA_ID_filtering {
            input:
                inputBAM = inputBAM,
                gtf_path = IDRefFree.mergedReffreeGTF,
                expr_file_path = Quant.mergedQuantExpr,
                Float threshold = 1.0
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/filtertranscripts:latest"
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
        File? refFreeUpdatedGTF = LRAA_ID_filtering.filtered_gtf
        File? refFreeQuant = Quant2.mergedQuantExprFile
        File? refFreeTracking = refFreeTracking.mergedQuantTracking
    }
}
