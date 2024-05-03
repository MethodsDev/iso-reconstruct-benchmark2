version 1.0

# This task uses IsoQuant version 3.3.0
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
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 256
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant@sha256:199395318713f631f29ef6f8b59debf3e0677a8a4f2590c7e9b7b941403d431f"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
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
            
            mv ~{OutDir}/ID_reffree/OUT/OUT.transcript_models.gtf ~{OutDir}/IsoQuant.gtf
            
            if [ -f ~{referenceAnnotation_reduced} ]; then
                /usr/local/src/IsoQuant/isoquant.py \
                --reference ~{referenceGenome} \
                --bam ~{inputBAM} \
                --genedb ~{referenceAnnotation_reduced} \
                --data_type ~{dataType} \
                --threads ~{numThreads} \
                --output ~{OutDir}/ID_reduced
                
            mv ~{OutDir}/ID_reduced/OUT/OUT.transcript_models.gtf ~{OutDir}/IsoQuant_reduced.gtf                              
            
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
            
            mv ~{OutDir}/Quant/OUT/OUT.transcript_counts.tsv ~{OutDir}/IsoQuant_quant.tsv 

        fi
    >>>

    command <<<
        if [ -f ~{inputBAM_with_polyA_for_IsoQuant} ] && [ -f ~{inputBAMIndex_with_polyA_for_IsoQuant} ]; then
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

            fi
        fi
    >>>

    output {
        File isoquantGTF = "~{OutDir}/IsoQuant.gtf"
        File? isoquantReducedGTF = "~{OutDir}/IsoQuant_reduced.gtf"
        File? isoquantCounts = "~{OutDir}/IsoQuant_quant.tsv"
        File? isoquantGTF_with_polyA = "~{OutDir}/IsoQuant_with_polyA.gtf"
        File? isoquantReducedGTF_with_polyA = "~{OutDir}/IsoQuant_reduced_with_polyA.gtf"
        File? isoquantCounts_with_polyA = "~{OutDir}/IsoQuant_quant_with_polyA.tsv"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
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
        Int cpu
        Int numThreads
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript
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
            ID_or_Quant_or_Both = ID_or_Quant_or_Both,
            cpu = cpu,
            numThreads = numThreads,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB,
            docker = docker,
            monitoringScript = monitoringScript
    }

    output {
        File isoquantGTF = isoquantTask.isoquantGTF
        File? isoquantReducedGTF = isoquantTask.isoquantReducedGTF
        File? isoquantCounts = isoquantTask.isoquantCounts
        File? isoquantGTF_with_polyA = isoquantTask.isoquantGTF_with_polyA
        File? isoquantReducedGTF_with_polyA = isoquantTask.isoquantReducedGTF_with_polyA
        File? isoquantCounts_with_polyA = isoquantTask.isoquantCounts_with_polyA
        File monitoringLog = isoquantTask.monitoringLog
    }
}