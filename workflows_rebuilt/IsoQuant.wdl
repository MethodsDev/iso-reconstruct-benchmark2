version 1.0

# This task uses IsoQuant version 3.4.0
task isoquantTask {
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
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:latest"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String OutDir = "IsoQuant_out"
    
    command <<<
        bash ~{monitoringScript} > monitoring.log &
        
        rm -rf ~{OutDir} && mkdir ~{OutDir} 

        if [ "~{ID_or_Quant_or_Both}" = "ID" ] || [ "~{ID_or_Quant_or_Both}" = "Both" ]; then
            /usr/local/src/IsoQuant/isoquant.py \
            --reference ~{referenceGenome} \
            --bam ~{inputBAM} \
            --data_type ~{dataType} \
            --threads ~{numThreads} \
            --output ~{OutDir}/ID_reffree
            
            mv ~{OutDir}/ID_reffree/OUT/OUT.transcript_models.gtf ~{OutDir}/isoquantGTF.gtf
            mv ~{OutDir}/ID_reffree/OUT/OUT.transcript_model_counts.tsv ~{OutDir}/isoquantGTFCounts.tsv




            if [ -f ~{referenceAnnotation_reduced} ]; then
                /usr/local/src/IsoQuant/isoquant.py \
                --reference ~{referenceGenome} \
                --bam ~{inputBAM} \
                --genedb ~{referenceAnnotation_reduced} \
                --data_type ~{dataType} \
                --threads ~{numThreads} \
                --output ~{OutDir}/ID_reduced
                
            mv ~{OutDir}/ID_reduced/OUT/OUT.transcript_models.gtf ~{OutDir}/isoquantReducedGTF.gtf                              
            mv ~{OutDir}/ID_reduced/OUT/OUT.transcript_model_counts.tsv ~{OutDir}/isoquantReducedGTFCounts.tsv
            fi
        fi

        if { [ "~{ID_or_Quant_or_Both}" = "Quant" ] || [ "~{ID_or_Quant_or_Both}" = "Both" ]; } && [ -f ~{referenceAnnotation_full} ]; then
            /usr/local/src/IsoQuant/isoquant.py \
            --reference ~{referenceGenome} \
            --bam ~{inputBAM} \
            --data_type ~{dataType} \
            --threads ~{numThreads} \
            --genedb ~{referenceAnnotation_full} \
            --output ~{OutDir}/Quant \
            --no_model_construction
            
            mv ~{OutDir}/Quant/OUT/OUT.transcript_counts.tsv ~{OutDir}/isoquantCounts.tsv
#            tar -czf IsoQuant_OUT.tar.gz ~{OutDir}/Quant

        fi

        if [ "~{inputBAM_with_polyA_for_IsoQuant}" != "" ] && [ "~{inputBAMIndex_with_polyA_for_IsoQuant}" != "" ]; then
            if [ "~{ID_or_Quant_or_Both}" = "ID" ] || [ "~{ID_or_Quant_or_Both}" = "Both" ]; then
                /usr/local/src/IsoQuant/isoquant.py \
                --reference ~{referenceGenome} \
                --bam ~{inputBAM_with_polyA_for_IsoQuant} \
                --data_type ~{dataType} \
                --threads ~{numThreads} \
                --output ~{OutDir}/ID_reffree_with_polyA
                
                mv ~{OutDir}/ID_reffree_with_polyA/OUT/OUT.transcript_models.gtf ~{OutDir}/IsoQuant_with_polyA.gtf

                if [ -f ~{referenceAnnotation_reduced} ]; then
                    /usr/local/src/IsoQuant/isoquant.py \
                    --reference ~{referenceGenome} \
                    --bam ~{inputBAM_with_polyA_for_IsoQuant} \
                    --genedb ~{referenceAnnotation_reduced} \
                    --data_type ~{dataType} \
                    --threads ~{numThreads} \
                    --output ~{OutDir}/ID_reduced_with_polyA
                    
                    mv ~{OutDir}/ID_reduced_with_polyA/OUT/OUT.transcript_models.gtf ~{OutDir}/IsoQuant_reduced_with_polyA.gtf

                fi
            fi

            if { [ "~{ID_or_Quant_or_Both}" = "Quant" ] || [ "~{ID_or_Quant_or_Both}" = "Both" ]; } && [ -f ~{referenceAnnotation_full} ]; then
                /usr/local/src/IsoQuant/isoquant.py \
                --reference ~{referenceGenome} \
                --bam ~{inputBAM_with_polyA_for_IsoQuant} \
                --data_type ~{dataType} \
                --threads ~{numThreads} \
                --genedb ~{referenceAnnotation_full} \
                --output ~{OutDir}/Quant_with_polyA \
                --no_model_construction
                
                mv ~{OutDir}/Quant_with_polyA/OUT/OUT.transcript_counts.tsv ~{OutDir}/IsoQuant_quant_with_polyA.tsv 
#                tar -czf IsoQuant_OUT_with_polyA.tar.gz ~{OutDir}/Quant_with_polyA

            fi
        fi
    >>>

    output {
        File? isoquantGTF = "~{OutDir}/isoquantGTF.gtf"
        File? isoquantReducedGTF = "~{OutDir}/isoquantReducedGTF.gtf"
        File? isoquantCounts = "~{OutDir}/isoquantCounts.tsv"

        File? isoquantGTF_with_polyA = "~{OutDir}/IsoQuant_with_polyA.gtf"
        File? isoquantReducedGTF_with_polyA = "~{OutDir}/IsoQuant_reduced_with_polyA.gtf"
        File? isoquantCounts_with_polyA = "~{OutDir}/IsoQuant_quant_with_polyA.tsv"
        File? isoquantCounts_OUT = "IsoQuant_OUT.tar.gz"
        File? isoquantCounts_with_polyA_OUT = "IsoQuant_OUT_with_polyA.tar.gz"

        File monitoringLog = "monitoring.log"
        File? isoquantReducedGTFCounts = "~{OutDir}/isoquantReducedGTFCounts.tsv"
        File? isoquantGTFCounts = "~{OutDir}/isoquantGTFCounts.tsv"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}


workflow isoquantWorkflow {
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
    }

    call isoquantTask {
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

    output {
        File? isoquantGTF = isoquantTask.isoquantGTF
        File? isoquantReducedGTF = isoquantTask.isoquantReducedGTF
        File? isoquantCounts = isoquantTask.isoquantCounts
        File? isoquantGTF_with_polyA = isoquantTask.isoquantGTF_with_polyA
        File? isoquantReducedGTF_with_polyA = isoquantTask.isoquantReducedGTF_with_polyA
        File? isoquantCounts_with_polyA = isoquantTask.isoquantCounts_with_polyA
        File? isoquantCounts_OUT = isoquantTask.isoquantCounts_OUT
        File? isoquantCounts_with_polyA_OUT = isoquantTask.isoquantCounts_with_polyA_OUT
        File monitoringLog = isoquantTask.monitoringLog
        File? isoquantReducedGTFCounts = isoquantTask.isoquantReducedGTFCounts
        File? isoquantGTFCounts = isoquantTask.isoquantGTFCounts
    }
}
