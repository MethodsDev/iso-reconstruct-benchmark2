version 1.0

import "Bambu.wdl" as bambuWorkflow
import "ESPRESSO.wdl" as espressoWorkflow
import "Flair.wdl" as flairWorkflow
import "FLAMES.wdl" as flamesWorkflow
import "IsoQuant.wdl" as isoquantWorkflow
import "IsoSeq.wdl" as isoseqWorkflow
import "Mandalorion.wdl" as mandalorionWorkflow
import "Mandalorion-fork.wdl" as mandalorionforkWorkflow
import "Oarfish.wdl" as oarfishWorkflow
import "Salmon.wdl" as salmonWorkflow
import "StringTie.wdl" as stringtieWorkflow
import "TALON.wdl" as talonWorkflow
import "LRAA.wdl" as lraaWorkflow
import "LRQuant.wdl" as lrquantWorkflow


task relocateOutputs {
    input {
        File? bambuReducedGTF
        File? bambuNDR1ReducedGTF
        File? bambuGTF
        File? bambuCounts
        File? espressoReducedGTF
        File? espressoCounts
        File? flairReducedGTF
        File? flairCounts
        File? flamesReducedGTF
        File? isoquantReducedGTF
        File? isoquantGTF
        File? isoquantCounts
        File? isoquantGTF_with_polyA
        File? isoquantReducedGTF_with_polyA
        File? isoquantCounts_with_polyA
        File? isoquantCounts_OUT
        File? isoquantCounts_with_polyA_OUT
        File? isoseqReducedGTF
        File? isoseqGTF
        File? mandalorionReducedGTF
        File? mandalorionGTF
        File? mandalorionforkReducedGTF
        File? mandalorionforkGTF
        File? oarfishCounts
        File? salmonCounts
        File? stringtieReducedGTF
        File? stringtieGTF
        File? stringtieCounts
        File? talonReducedGTF
        File? lraaGTF
        File? lraaReducedGTF
        File? lraaCounts
        File? lraaCounts_noEM
        File? lraa_quant_tracking
        File? lraa_quant_tracking_noEM
        File? gffcompareCounts
        File? lrquantCounts
        File? lrquantOUT

        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/espresso@sha256:f538303f6457c55e7b3c2a45081e6d8e3053e6f76e56bc65631b7f4aa290b026"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
        Int cpu = 16
        Int memoryGB = 256
        Int diskSizeGB = 500
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        mkdir ID_reduced ID Quant All_Outputs_Relocated
    
        # Define arrays of files for each directory
        reduced_files=("~{bambuReducedGTF}" "~{bambuNDR1ReducedGTF}" "~{espressoReducedGTF}" "~{flairReducedGTF}" "~{flamesReducedGTF}" "~{isoquantReducedGTF}" "~{isoquantReducedGTF_with_polyA}" "~{isoseqReducedGTF}" "~{mandalorionReducedGTF}" "~{mandalorionforkReducedGTF}" "~{stringtieReducedGTF}" "~{talonReducedGTF}" "~{lraaReducedGTF}")
        id_files=("~{bambuGTF}" "~{isoquantGTF}" "~{isoquantGTF_with_polyA}" "~{isoseqGTF}" "~{mandalorionGTF}" "~{mandalorionforkGTF}" "~{stringtieGTF}" "~{lraaGTF}")
        quant_files=("~{bambuCounts}" "~{espressoCounts}" "~{flairCounts}" "~{isoquantCounts}" "~{isoquantCounts_with_polyA}" "~{oarfishCounts}" "~{salmonCounts}" "~{stringtieCounts}" "~{lraaCounts}" "~{lraaCounts_noEM}" "~{lraa_quant_tracking}" "~{lraa_quant_tracking_noEM}" "~{gffcompareCounts}" "~{lrquantCounts}")
    
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

        mv ID_reduced ID Quant All_Outputs_Relocated/
        tar -czf All_Outputs_Relocated.tar.gz All_Outputs_Relocated/

#        if [[ -n "~{lrquantOUT}" ]]; then
#            mv ~{lrquantOUT} LRQuant_OUT.tar.gz
#        fi
    
#        if [[ -n "~{isoquantCounts_OUT}" ]]; then
#            mv ~{isoquantCounts_OUT} isoquantCounts_OUT.tar.gz
#        fi
    
#        if [[ -n "~{isoquantCounts_with_polyA_OUT}" ]]; then
#            mv ~{isoquantCounts_with_polyA_OUT} isoquantCounts_with_polyA_OUT.tar.gz
    >>>

    output {
        File? relocated_files = "All_Outputs_Relocated.tar.gz"
#        File? lrquant_files = "LRQuant_OUT.tar.gz"
#        File? isoquantCounts_OUT = "isoquantCounts_OUT.tar.gz"
#        File? isoquantCounts_with_polyA_OUT = "isoquantCounts_with_polyA_OUT.tar.gz"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
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
        Boolean? LRAA_no_norm = false
        Boolean runBambu = true
        Boolean runEspresso = true
        Boolean runFlair = true
        Boolean runFlames = true
        Boolean runIsoquant = true
        Boolean runIsoseq = true
        Boolean runMandalorionfork = true
        Boolean runMandalorion = true
        Boolean runOarfish = true
        Boolean runSalmon = true
        Boolean runStringtie = true
        Boolean runTalon = true
        Boolean runLraa = true
        Boolean runLrquant = true
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
    call isoseqWorkflow.isoseqWorkflow as isoseq {
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

if (runMandalorionfork) {
    call mandalorionforkWorkflow.mandalorionforkWorkflow as mandalorionfork {
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
    

if (runLraa) {
    call lraaWorkflow.lraaWorkflow as lraa {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both,
            LRAA_no_norm = LRAA_no_norm
    }
}

if (runLrquant) {
    call lrquantWorkflow.lrquantWorkflow as lrquant {
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
        isoquantReducedGTF = isoquant.isoquantReducedGTF,
        isoquantGTF = isoquant.isoquantGTF,
        isoquantCounts = isoquant.isoquantCounts,
        isoquantGTF_with_polyA = isoquant.isoquantGTF_with_polyA,
        isoquantReducedGTF_with_polyA = isoquant.isoquantReducedGTF_with_polyA,
        isoquantCounts_with_polyA = isoquant.isoquantCounts_with_polyA,
 #       isoquantCounts_OUT = isoquant.isoquantCounts_OUT,
 #       isoquantCounts_with_polyA_OUT = isoquant.isoquantCounts_with_polyA_OUT,  
        isoseqReducedGTF = isoseq.isoseqReducedGTF,
        isoseqGTF = isoseq.isoseqGTF,
        mandalorionReducedGTF = mandalorion.mandalorionReducedGTF,
        mandalorionGTF = mandalorion.mandalorionGTF,
        mandalorionforkReducedGTF = mandalorionfork.mandalorionforkReducedGTF,
        mandalorionforkGTF = mandalorionfork.mandalorionforkGTF,
        oarfishCounts = oarfish.oarfishCounts,
        salmonCounts = salmon.salmonCounts,
        stringtieReducedGTF = stringtie.stringtieReducedGTF,
        stringtieGTF = stringtie.stringtieGTF,
        stringtieCounts = stringtie.stringtieCounts,
        talonReducedGTF = talon.talonReducedGTF,
        lraaGTF = lraa.lraaGTF,
        lraaReducedGTF = lraa.lraaReducedGTF,
        lraaCounts = lraa.lraaCounts,
        lraaCounts_noEM = lraa.lraaCounts_noEM,
        lraa_quant_tracking = lraa.lraa_quant_tracking,
        lraa_quant_tracking_noEM = lraa.lraa_quant_tracking_noEM
#        gffcompareCounts = lrquant.gffcompareCounts,
#        lrquantCounts = lrquant.lrquantCounts,
#        lrquantOUT = lrquant.lrquantOUT
}
}


