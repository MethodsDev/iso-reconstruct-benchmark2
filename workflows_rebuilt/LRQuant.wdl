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
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 250
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/lrquant:latest"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
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
            -p LRQuant_OUT \
            -t ~{numThreads}

        fi
    >>>

    output {
        File? gffcompareCounts = "Gffcompare_quant.tsv"
        File? lrquantCounts = "LRQuant_OUT_expression_matrix.tsv"
        File? lrquantOUT = "LRQuant_OUT.tar.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
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
        File? lrquantCounts =lrquantTask.lrquantCounts
        File? gffcompareCounts = lrquantTask.gffcompareCounts
        File? lrquantOUT = lrquantTask.lrquantOUT       
        File monitoringLog = lrquantTask.monitoringLog
    }
}
