version 1.0

# This task uses Isosceles version
task isoscelesTask {
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
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isosceles:latest"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }
    
    String OutDir = "Isosceles_out"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        if [[ "~{ID_or_Quant_or_Both}" == "ID" ]]; then
            if [[ -n "~{referenceAnnotation_reduced}" ]]; then
                isosceles -b ~{inputBAM} \
                -i ~{referenceAnnotation_reduced} \
                -f ~{referenceGenome} -n ~{numThreads} -m ~{dataType} -t ID
            fi
        fi

        if [[ "~{ID_or_Quant_or_Both}" == "Quant" ]]; then
            if [[ -n "~{referenceAnnotation_full}" ]]; then
                isosceles -b ~{inputBAM} \
                -q ~{referenceAnnotation_full} \
                -f ~{referenceGenome} -n ~{numThreads} -m ~{dataType} -t Quant
            fi
        fi

        if [[ "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            if [[ -n "~{referenceAnnotation_reduced}" && -n "~{referenceAnnotation_full}" ]]; then
                isosceles -b ~{inputBAM} \
                -i ~{referenceAnnotation_reduced} \
                -q ~{referenceAnnotation_full} \
                -f ~{referenceGenome} -n ~{numThreads} -m ~{dataType} -t Both
            fi
        fi
    >>>
    
    output {
        File? isoscelesReducedGTF = "Isosceles_de_novo_loose.gtf"
        File? isoscelesStrictReducedGTF = "Isosceles_de_novo_strict.gtf"
        File? isoscelesCounts = "Isosceles_quant.txt"
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

workflow isoscelesWorkflow {
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

    call isoscelesTask {
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
        File? isoscelesReducedGTF = isoscelesTask.isoscelesReducedGTF
        File? isoscelesStrictReducedGTF = isoscelesTask.isoscelesStrictReducedGTF
        File? isoscelesCounts = isoscelesTask.isoscelesCounts
        File monitoringLog = isoscelesTask.monitoringLog
    }
}
