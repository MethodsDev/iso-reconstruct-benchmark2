version 1.0

import "ID_ref_free.wdl" as IDRefFree
import "ID_ref_guided.wdl" as IDRefGuided
import "Quant.wdl" as Quant
import "Filtering.wdl" as Filtering

workflow CombinedWorkflow {
    input {
        File inputBAM
        File referenceGenome
        File referenceGTF
        String mode
        # Additional inputs as required by the individual tasks
    }

    # Conditional execution based on the mode
    if (mode == "ID_ref_free_Quant_mode") {
        call IDRefFree.IDRefFreeTask {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome
                # Add other necessary inputs for IDRefFreeTask
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
