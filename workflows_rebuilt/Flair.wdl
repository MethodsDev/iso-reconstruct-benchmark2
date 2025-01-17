version 1.0

# This task uses Flair version 2.0.0
task flairTask {
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
        Int diskSizeGB = 512
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/flair:latest"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String flairPrefix = "Flair"
    String OutDir = "Flair_out"

    command <<<
        bash ~{monitoringScript} > monitoring.log &
        
        mkdir -p ~{OutDir}
        cd ~{OutDir}
        
        samtools bam2fq ~{inputBAM} > temp.fastq
        bam2Bed12 -i ~{inputBAM} > "~{flairPrefix}.bed"

        if [ "~{ID_or_Quant_or_Both}" = "ID" -o "~{ID_or_Quant_or_Both}" = "Both" ]; then
            if [[ -n "~{referenceAnnotation_reduced}" ]]; then
                flair correct \
                -q "~{flairPrefix}.bed" \
                --genome ~{referenceGenome} \
                --gtf ~{referenceAnnotation_reduced} \
                -o ~{flairPrefix} \
                -t ~{numThreads}

                flair collapse \
                -g ~{referenceGenome} \
                -f ~{referenceAnnotation_reduced} \
                -r temp.fastq \
                -q "~{flairPrefix}_all_corrected.bed" \
                -o ~{flairPrefix} \
                -t ~{numThreads} #\
       #         --stringent --check_splice --generate_map --annotation_reliant generate

                mv Flair.isoforms.gtf flairReducedGTF.gtf


                gffread flairReducedGTF.gtf -g ~{referenceGenome} -w transcriptome2.fa
                sample2=("flair" "condition1" "batch1" "temp.fastq")
                manifest_filename2="flair_manifest2.tsv"
    
                for i in "${sample2[@]}"; do
                    echo -e "$i\t\c" >> $manifest_filename2
                done
                echo "" >> $manifest_filename2

                flair quantify -r $manifest_filename2 -i transcriptome2.fa -o ReducedQuant -t ~{numThreads}
                mv ReducedQuant.counts.tsv flairReducedGTFCounts.tsv
            fi
        fi

        if [ "~{ID_or_Quant_or_Both}" = "Quant" -o "~{ID_or_Quant_or_Both}" = "Both" ]; then
        
            gffread ~{referenceAnnotation_full} -g ~{referenceGenome} -w transcriptome.fa
            sample=("flair" "condition1" "batch1" "temp.fastq")
            manifest_filename="flair_manifest.tsv"

            for i in "${sample[@]}"; do
                echo -e "$i\t\c" >> $manifest_filename
            done
            echo "" >> $manifest_filename
            
            flair quantify -r $manifest_filename -i transcriptome.fa -o Quant -t ~{numThreads}
            
            mv Quant.counts.tsv flairCounts.tsv
        fi
    >>>

    output {
        File? flairReducedGTF = "~{OutDir}/flairReducedGTF.gtf"
        File? flairCounts = "~{OutDir}/flairCounts.tsv"
        File monitoringLog = "monitoring.log"
        File? flairReducedGTFCounts = "~{OutDir}/flairReducedGTFCounts.tsv"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}


workflow flairWorkflow {
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

    call flairTask {
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
        File? flairReducedGTF = flairTask.flairReducedGTF
        File? flairCounts = flairTask.flairCounts
        File monitoringLog = flairTask.monitoringLog
    }
}
