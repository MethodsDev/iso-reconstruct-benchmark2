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
}


    Int diskSizeGB = 1024
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    String OutDir = "LRAA_out"
    String docker_filtering


    if (mode == "ID_ref_free_Quant_mode") {

        call IDRefFree.lraaWorkflow {
            input:
                File inputBAM
                File referenceGenome
                Int numThreads = numThreads
                Int memoryGB = memoryGB
                Int diskSizeGB = diskSizeGB
                String docker = docker
                String main_chromosomes = main_chromosomes
                Boolean? LRAA_no_norm = LRAA_no_norm
        }





        call Quant.QuantTask {
            input:
                inputBAM = inputBAM,
                referenceGTF = IDRefFree.IDRefFreeTask.updatedGTF
                # Add other necessary inputs for QuantTask
        }

        call Filtering.FilteringTask {
            input:
                inputGTF = IDRefFree.IDRefFreeTask.updatedGTF,
                quantResults = Quant.QuantTask.quantResults
                # Add other necessary inputs for FilteringTask
        }

        call Quant.QuantTask as QuantAfterFiltering {
            input:
                inputBAM = inputBAM,
                referenceGTF = Filtering.FilteringTask.filteredGTF
                # Add other necessary inputs for QuantTask
        }
    }

    if (mode == "ID_ref_guided_Quant_mode") {
        call IDRefGuided.IDRefGuidedTask {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome
                # Add other necessary inputs for IDRefGuidedTask
        }

        call Quant.QuantTask as QuantForRefGuided {
            input:
                inputBAM = inputBAM,
                referenceGTF = IDRefGuided.IDRefGuidedTask.updatedGTF
                # Add other necessary inputs for QuantTask
        }

        call Filtering.FilteringTask as FilteringForRefGuided {
            input:
                inputGTF = IDRefGuided.IDRefGuidedTask.updatedGTF,
                quantResults = QuantForRefGuided.quantResults
                # Add other necessary inputs for FilteringTask
        }

        call Quant.QuantTask as QuantAfterFilteringForRefGuided {
            input:
                inputBAM = inputBAM,
                referenceGTF = FilteringForRefGuided.filteredGTF
                # Add other necessary inputs for QuantTask
        }
    }

    if (mode == "Quant_only") {
        call Quant.QuantTask as QuantOnly {
            input:
                inputBAM = inputBAM,
                referenceGTF = referenceGTF
                # Add other necessary inputs for QuantTask
        }
    }

    # Define workflow outputs
    output {
        File? refFreeUpdatedGTF = IDRefFree.IDRefFreeTask.updatedGTF
        File? refGuidedUpdatedGTF = IDRefGuided.IDRefGuidedTask.updatedGTF
        File? quantResultsFirstPass = select_first([Quant.QuantTask.quantResults, QuantForRefGuided.quantResults])
        File? filteredGTF = select_first([Filtering.FilteringTask.filteredGTF, FilteringForRefGuided.filteredGTF])
        File? finalQuantResults = select_first([QuantAfterFiltering.quantResults, QuantAfterFilteringForRefGuided.quantResults, QuantOnly.quantResults])
    }
}
