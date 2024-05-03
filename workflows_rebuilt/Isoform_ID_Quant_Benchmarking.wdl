version 1.0

import "Bambu.wdl" as bambuWorkflow
import "ESPRESSO.wdl" as espressoWorkflow
import "Flair.wdl" as flairWorkflow
import "FLAMES-py.wdl" as flamesWorkflow
import "IsoQuant.wdl" as isoquantWorkflow
import "IsoSeq.wdl" as isoseqWorkflow
import "Mandalorion.wdl" as mandalorionWorkflow
import "Oarfish.wdl" as oarfishWorkflow
import "Salmon.wdl" as salmonWorkflow
import "StringTie.wdl" as stringtieWorkflow
import "TALON.wdl" as talonWorkflow

task relocateOutputs {
    input {
        ?File bambuReducedGTF
        ?File bambuNDR1ReducedGTF
        ?File bambuGTF
        ?File bambuCounts
        ?File espressoReducedGTF
        ?File espressoCounts
        ?File flairReducedGTF
        ?File flairCounts
        ?File flamesReducedGTF
        ?File isoquantReducedGTF
        ?File isoquantGTF
        ?File isoquantCounts
        ?File isoquantGTF_with_polyA
        ?File isoquantReducedGTF_with_polyA
        ?File isoquantCounts_with_polyA
        ?File isoseqReducedGTF
        ?File isoseqGTF
        ?File mandalorionReducedGTF
        ?File mandalorionGTF
        ?File oarfishCounts
        ?File salmonCounts
        ?File stringtieReducedGTF
        ?File stringtieGTF
        ?File stringtieCounts
        ?File talonReducedGTF
    }

    command {
        mkdir ID_reduced ID Quant

        # Define arrays of files for each directory
        reduced_files=(~{bambuReducedGTF} ~{bambuNDR1ReducedGTF} ~{espressoReducedGTF} ~{flairReducedGTF} ~{flamesReducedGTF} ~{isoquantReducedGTF} ~{isoquantReducedGTF_with_polyA} ~{isoseqReducedGTF} ~{mandalorionReducedGTF} ~{stringtieReducedGTF} ~{talonReducedGTF})
        id_files=(~{bambuGTF} ~{isoquantGTF} ~{isoquantGTF_with_polyA} ~{isoseqGTF} ~{mandalorionGTF} ~{stringtieGTF})
        quant_files=(~{bambuCounts} ~{espressoCounts} ~{flairCounts} ~{isoquantCounts} ~{isoquantCounts_with_polyA} ~{oarfishCounts} ~{salmonCounts} ~{stringtieCounts})

        # Loop over the files for each directory
        for file in "${reduced_files[@]}"; do
          [ -f "$file" ] && mv "$file" ID_reduced/
        done

        for file in "${id_files[@]}"; do
          [ -f "$file" ] && mv "$file" ID/
        done

        for file in "${quant_files[@]}"; do
          [ -f "$file" ] && mv "$file" Quant/
        done
    }

    output {
        Directory ID_reduced = "ID_reduced"
        Directory ID = "ID"
        Directory Quant = "Quant"
    }

    runtime {
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 10 HDD"
    }
}

workflow LongReadRNABenchmark {
    input {
        File inputBAM
        File inputBAMIndex
        File? inputBAM_with_polyA_for_IsoQuant
        File? inputBAMIndex_with_polyA_for_IsoQuant
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String dataType
        String ID_or_Quant_or_Both
        Boolean runBambu = true
        Boolean runEspresso = true
        Boolean runFlair = true
        Boolean runFlames = true
        Boolean runIsoquant = true
        Boolean runIsoseq = true
        Boolean runMandalorion = true
        Boolean runOarfish = true
        Boolean runSalmon = true
        Boolean runStringtie = true
        Boolean runTalon = true
    }
if (runBambu) {
    call bambuWorkflow.bambuWorkflow as bambu {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }
}

if (runEspresso) {
    call espressoWorkflow.espressoWorkflow as espresso {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }
}

if (runFlair) {
    call flairWorkflow.flairWorkflow as flair {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }
}

if (runFlames) {
    call flamesWorkflow.flamesWorkflow as flames {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }
}

if (runIsoquant) {
    call isoquantWorkflow.isoquantWorkflow as isoquant {
        input:   
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            inputBAM_with_polyA_for_IsoQuant = inputBAM_with_polyA_for_IsoQuant,
            inputBAMIndex_with_polyA_for_IsoQuant = inputBAMIndex_with_polyA_for_IsoQuant,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }
}

if (runIsoseq) {
    call isoseqWorkflow.isoseqWorkflow as isoeq {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }
}

if (runMandalorion) {
    call mandalorionWorkflow.mandalorionWorkflow as mandalorion {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }
}

if (runOarfish) {
    call oarfishWorkflow.oarfishWorkflow as oarfish {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }
}

if (runSalmon) {
    call salmonWorkflow.salmonWorkflow as salmon {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }
}

if (runStringtie) {
    call stringtieWorkflow.stringtieWorkflow as stringtie {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }
}

if (runTalon) {
    call talonWorkflow.talonWorkflow as talon {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both
    }
}
    
    
call relocateOutputs {
    input:
        bambuReducedGTF = bambu.bambuReducedGTF,
        bambuNDR1ReducedGTF = bambu.bambuNDR1ReducedGTF,
        bambuGTF = bambu.bambuGTF,
        bambuCounts = bambu.bambuCounts,
        espressoReducedGTF = espresso.espressoReducedGTF,
        espressoCounts = espresso.espressoCounts,
        flairReducedGTF = flair.flairReducedGTF,
        flairCounts = flair.flairCounts,
        flamesReducedGTF = flames.flamesReducedGTF,
        flamesCounts = flames.flamesCounts,
        isoquantReducedGTF = isoquant.isoquantReducedGTF,
        isoquantGTF = isoquant.isoquantGTF,
        isoquantCounts = isoquant.isoquantCounts,
        isoquantGTF_with_polyA = isoquant.isoquantGTF_with_polyA,
        isoquantReducedGTF_with_polyA = isoquant.isoquantReducedGTF_with_polyA,
        isoquantCounts_with_polyA = isoquant.isoquantCounts_with_polyA,
        isoseqReducedGTF = isoseq.isoseqReducedGTF,
        isoseqGTF = isoseq.isoseqGTF,
        mandalorionReducedGTF = mandalorion.mandalorionReducedGTF,
        mandalorionGTF = mandalorion.mandalorionGTF,
        oarfishCounts = oarfish.oarfishCounts,
        salmonCounts = salmon.salmonCounts,
        stringtieReducedGTF = stringtie.stringtieReducedGTF,
        stringtieGTF = stringtie.stringtieGTF,
        stringtieCounts = stringtie.stringtieCounts,
        talonReducedGTF = talon.talonReducedGTF
}
}


