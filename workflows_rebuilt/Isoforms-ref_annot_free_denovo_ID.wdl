version 1.0

import "Bambu.wdl" as bambuWorkflow
import "Flair.wdl" as flairWorkflow
import "FLAMES.wdl" as flamesWorkflow
import "IsoQuant.wdl" as isoquantWorkflow
import "IsoSeq.wdl" as isoseqWorkflow
import "Mandalorion.wdl" as mandalorionWorkflow
import "StringTie.wdl" as stringtieWorkflow
import "TALON.wdl" as talonWorkflow
import "https://raw.githubusercontent.com/MethodsDev/LongReadAlignmentAssembler/refs/heads/main/WDL/LRAA.wdl" as lraaWorkflow
import "Isosceles.wdl" as isoscelesWorkflow


workflow QuantOnly_wf {

    input {
        String sample_id
        File inputBAM
        File inputBAMIndex
        File referenceGenomeFasta
        File referenceGenomeIndex

        String data_type
        String oarfish_seq_tech
        Boolean LRAA_LowFi
        String? LRAA_main_chromosomes
        
        Boolean runBambu = true
        Boolean runIsoquant = true
        Boolean runIsoseq = true
        Boolean runMandalorion = true
        Boolean runStringtie = true
        Boolean runLRAA = true

        Boolean AggregateOutputs = true
        
    }
    
    if (runBambu) {
        call bambuWorkflow.bambuWorkflow as bambu {
            input:
                sample_id = sample_id,
                inputBAM = inputBAM,
                inputBAMIndex = inputBAMIndex,
                referenceGenomeFasta = referenceGenomeFasta,
                referenceGenomeIndex = referenceGenomeIndex,
                quant_only = false
        }
    }

    if (runIsoquant) {
        call isoquantWorkflow.isoquantWorkflow as isoquant {
            input:
                sample_id = sample_id,
                inputBAM = inputBAM,
                inputBAMIndex = inputBAMIndex,
                referenceGenomeFasta = referenceGenomeFasta,
                referenceGenomeIndex = referenceGenomeIndex,
                data_type = data_type,
                quant_only = false
        }
    }

    if (runIsoseq) {
        call isoseqWorkflow.isoseqWorkflow as isoseq {
            input:
                sample_id = sample_id,
                inputBAM = inputBAM,
                inputBAMIndex = inputBAMIndex,
                referenceGenomeFasta = referenceGenomeFasta,
                referenceGenomeIndex = referenceGenomeIndex
        }
    }

    if (runMandalorion) {
        call mandalorionWorkflow.mandalorionWorkflow as mandalorion {
            input:
                sample_id = sample_id,
                inputBAM = inputBAM,
                inputBAMIndex = inputBAMIndex,
                referenceGenomeFasta = referenceGenomeFasta,
                referenceGenomeIndex = referenceGenomeIndex
        }
    }

    if (runStringtie) {
        call stringtieWorkflow.stringtieWorkflow as stringtie {
            input:
                sample_id = sample_id,
                inputBAM = inputBAM,
                inputBAMIndex = inputBAMIndex,
                referenceGenomeFasta = referenceGenomeFasta,
                referenceGenomeIndex = referenceGenomeIndex,
                quant_only = false
        }
    }

    if (runLRAA) {
        call lraaWorkflow.LRAA_wf as lraa {
            input:
                sample_id = sample_id,
                inputBAM = inputBAM,
                referenceGenome = referenceGenomeFasta,
                LowFi = LRAA_LowFi,
                main_chromosomes = LRAA_main_chromosomes,
                quant_only = false
        }
    }


    if (AggregateOutputs) {
        call AggregateOutputs {
            input:

            sample_id = sample_id,

            bambu_quant = bambu.bambu_quant,
            bambu_gtf = bambu.bambu_gtf,

            isoquant_gtf = isoquant.isoquant_gtf,
            isoquant_counts = isoquant.isoquant_counts,

            isoseq_gtf = isoseq.isoseq_ref_filtered_GTF,
            isoseq_counts = isoseq.isoseq_counts,

            mandalorion_gtf = mandalorion.mandalorion_gtf,
            mandalorion_quant = mandalorion.mandalorion_quant,

            stringtie_gtf = stringtie.stringtie_gtf,
            stringtie_quant = stringtie.stringtie_quant,

            LRAA_gtf = lraa.mergedGTF,
            LRAA_quant = lraa.mergedQuantExpr,

        }
    }
}



task AggregateOutputs {
    input {
        
        String sample_id
        
        File? bambu_quant
        File? bambu_gtf

        File? isoquant_gtf
        File? isoquant_counts

        File? isoseq_gtf
        File? isoseq_counts

        File? mandalorion_gtf
        File? mandalorion_quant

        File? stringtie_gtf
        File? stringtie_quant

        File? LRAA_gtf
        File? LRAA_quant
        
        String docker = "ubuntu"

        Int cpu = 4
        Int memoryGB = 16
        Int diskSizeGB = 500
    }

    
    String output_dir = "All_RefAnnotFreeDenovoID_Outputs_Aggregated-~{sample_id}"

    
    command <<<

        set -ex

         find . -type f | grep "~{sample_id}" | xargs -n1 -I__ mv __ . 
        
        ls -ltr
        
        mkdir ~{output_dir}

        mv ~{sample_id}* ~{output_dir}/
        
        tar -czvf ~{output_dir}.tar.gz ~{output_dir}/

        ls -ltr

        echo Done.
   

    >>>

    output {
        File? aggregated_files = "~{output_dir}.tar.gz"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}
