version 1.0

# This task uses LRQuant version 0.2
task lrquantTask {
    input {
        File inputBAM
        File inputBAMIndex
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
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/lrquant:latest"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
    }
    String OutDir = "LRQuant_out"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        if [[ "~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            if [[ -z "~{referenceAnnotation_full}" ]]; then
                echo "Full annotation is not provided. Please provide the full annotation."
                exit 1
            fi  
            mkdir -p ~{OutDir}


            samtools bam2fq --threads ~{numThreads} ~{inputBAM} > LRQuant_tmp.fq
            LRQuant -r LRQuant_tmp.fq \
            -g ~{referenceGenome} \
            -a ~{referenceAnnotation_full} \
            -p LRQuant_OUT


            mv LRQuant_OUT/gffcompare_out/LRQuant_OUT.gffcompare.tsv Gffcompare_quant.tsv
            mv LRQuant_OUT/lrquant_out/LRQuant_OUT.lrquant.tsv LRQuant_quant.tsv
            tar -czf LRQuant_OUT.tar.gz LRQuant_OUT/
        fi
    >>>

    output {
        File? gffcompareCounts = "LRQuant_quant.tsv"
        File? lrquantCounts = "Gffcompare_quant.tsv"
        File? lrquantOUT = "LRQuant_OUT.zip"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}


workflow lrquantWorkflow {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String dataType
        String ID_or_Quant_or_Both
    }

    call lrquantTask {
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

    output {
        File? gffcompareCounts = "LRQuant_quant.tsv"
        File? lrquantCounts = "Gffcompare_quant.tsv"
        File? lrquantOUT = "LRQuant_OUT.zip"        
        File monitoringLog = lrquantTask.monitoringLog
    }
}
