version 1.0

# This task uses Bambu version 3.4.0
task bambuTask {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String dataType
        String ID_or_Quant_or_Both
        String Reffree_or_Refguided_or_Both
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 250
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/bambu@sha256:04601e3dd0d7f9c2008c5f430da7aef62277dc8ffd50ecb356d4b48ba962e6c7"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String OutDir = "Bambu_out"

    command <<<

    bash ~{monitoringScript} > monitoring.log &
    mkdir ~{OutDir}
    mkdir ~{OutDir}/ID_reduced
    mkdir ~{OutDir}/ID_ndr1_reduced
    mkdir ~{OutDir}/ID_reffree
    mkdir ~{OutDir}/Quant

    if [[ "~{ID_or_Quant_or_Both}" != "ID" && -z "~{referenceAnnotation_full}" ]]; then
        echo "Error: referenceAnnotation_full must be provided if ID_or_Quant_or_Both is not equal to ID."
        exit 1
    fi

    if [ "~{ID_or_Quant_or_Both}" == "ID" ] || [ "~{ID_or_Quant_or_Both}" == "Both" ]; then
        if [ "~{Reffree_or_Refguided_or_Both}" == "Reffree" ] || [ "~{Reffree_or_Refguided_or_Both}" == "Both" ]; then
            Rscript -<< EOF
            library(bambu)
            fa.file <- "~{referenceGenome}"
            lr.bam <- "~{inputBAM}"
            lr.se <- bambu(reads = lr.bam, rcOutDir = '~{OutDir}/ID_reffree', annotations = NULL, genome = fa.file, quant = FALSE, NDR = 1, ncore = ~{numThreads})
            writeToGTF(lr.se, "~{OutDir}/ID_reffree/bambuGTF.gtf")
EOF
            
            mv ~{OutDir}/ID_reffree/bambuGTF.gtf  ~{OutDir}/bambuGTF.gtf

            Rscript -<< EOF
            library(bambu)
            test.bam <- "~{inputBAM}"
            fa.file <- "~{referenceGenome}"
            gtf.file <- "~{OutDir}/bambuGTF.gtf"
            se.quantOnly <- bambu(reads = test.bam, annotations = gtf.file, genome = fa.file, discovery = FALSE)
            writeBambuOutput(se.quantOnly, path = "~{OutDir}/Quant_ReffreeGTF")
EOF
            mv ~{OutDir}/Quant_ReffreeGTF/counts_transcript.txt ~{OutDir}/bambuGTFCounts.txt
        fi

        if [ "~{Reffree_or_Refguided_or_Both}" == "Refguided" ] || [ "~{Reffree_or_Refguided_or_Both}" == "Both" ]; then
            if [ "~{referenceAnnotation_reduced}" != "" ]; then
                Rscript -<< EOF
                library(bambu)
                fa.file <- "~{referenceGenome}"
                gtf.file <- "~{referenceAnnotation_reduced}"
                bambuAnnotations <- prepareAnnotations(gtf.file)
                lr.bam <- "~{inputBAM}"
                lr.se <- bambu(reads = lr.bam, rcOutDir = "~{OutDir}/ID_reduced", annotations = bambuAnnotations, genome = fa.file, ncore = ~{numThreads})
                writeBambuOutput(lr.se, path = "~{OutDir}/ID_reduced")
EOF

                awk ' $3 >= 1 ' ~{OutDir}/ID_reduced/counts_transcript.txt | sort -k3,3n > ~{OutDir}/ID_reduced/expressed_annotations.gtf.counts
                cut -f1 ~{OutDir}/ID_reduced/expressed_annotations.gtf.counts > ~{OutDir}/ID_reduced/expressed_transcripts.txt
                grep -Ff ~{OutDir}/ID_reduced/expressed_transcripts.txt ~{OutDir}/ID_reduced/extended_annotations.gtf > ~{OutDir}/bambuReducedGTF.gtf
                mv ~{OutDir}/ID_reduced/expressed_annotations.gtf.counts ~{OutDir}/bambuReducedGTFCounts.txt

                Rscript -<< EOF
                library(bambu)
                fa.file <- "~{referenceGenome}"
                gtf.file <- "~{referenceAnnotation_reduced}"
                bambuAnnotations <- prepareAnnotations(gtf.file)
                lr.bam <- "~{inputBAM}"
                lr.se <- bambu(reads = lr.bam, rcOutDir = "~{OutDir}/ID_ndr1_reduced", annotations = bambuAnnotations, genome = fa.file, ncore = ~{numThreads}, NDR = 1)
                writeBambuOutput(lr.se, path = "~{OutDir}/ID_ndr1_reduced")
EOF

                awk ' $3 >= 1 ' ~{OutDir}/ID_ndr1_reduced/counts_transcript.txt | sort -k3,3n > ~{OutDir}/ID_ndr1_reduced/expressed_annotations.gtf.counts
                cut -f1 ~{OutDir}/ID_ndr1_reduced/expressed_annotations.gtf.counts > ~{OutDir}/ID_ndr1_reduced/expressed_transcripts.txt
                grep -Ff ~{OutDir}/ID_ndr1_reduced/expressed_transcripts.txt ~{OutDir}/ID_ndr1_reduced/extended_annotations.gtf > ~{OutDir}/bambuNDR1ReducedGTF.gtf
                mv ~{OutDir}/ID_ndr1_reduced/expressed_annotations.gtf.counts ~{OutDir}/bambuNDR1ReducedGTFCounts.txt
            fi
        fi
    fi

    if [ "~{referenceAnnotation_full}" != "" ] && ([ "~{ID_or_Quant_or_Both}" == "Quant" ] || [ "~{ID_or_Quant_or_Both}" == "Both" ]); then
        Rscript -<< EOF
        library(bambu)
        test.bam <- "~{inputBAM}"
        fa.file <- "~{referenceGenome}"
        gtf.file <- "~{referenceAnnotation_full}"
        se.quantOnly <- bambu(reads = test.bam, annotations = gtf.file, genome = fa.file, discovery = FALSE)
        writeBambuOutput(se.quantOnly, path = "~{OutDir}/Quant")
EOF

        mv ~{OutDir}/Quant/counts_transcript.txt ~{OutDir}/bambuCounts.txt
    fi

    >>>

    output {
        File? bambuReducedGTF = "~{OutDir}/bambuReducedGTF.gtf"
        File? bambuNDR1ReducedGTF = "~{OutDir}/bambuNDR1ReducedGTF.gtf"
        File? bambuGTF = "~{OutDir}/bambuGTF.gtf"
        File? bambuCounts = "~{OutDir}/bambuCounts.txt"
        File monitoringLog = "monitoring.log"
        File? bambuReducedGTFCounts = "~{OutDir}/bambuReducedGTFCounts.txt"
        File? bambuNDR1ReducedGTFCounts = "~{OutDir}/bambuNDR1ReducedGTFCounts.txt"
        File? bambuGTFCounts = "~{OutDir}/bambuGTFCounts.txt"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
    }
}


workflow bambuWorkflow {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String dataType
        String ID_or_Quant_or_Both
        String Reffree_or_Refguided_or_Both
    }

    call bambuTask {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both,
            Reffree_or_Refguided_or_Both = Reffree_or_Refguided_or_Both
    }

    output {
        File? bambuReducedGTF = bambuTask.bambuReducedGTF
        File? bambuNDR1ReducedGTF = bambuTask.bambuNDR1ReducedGTF
        File? bambuGTF = bambuTask.bambuGTF
        File? bambuCounts = bambuTask.bambuCounts
        File monitoringLog = bambuTask.monitoringLog
        File? bambuReducedGTFCounts = bambuTask.bambuReducedGTFCounts
        File? bambuNDR1ReducedGTFCounts = bambuTask.bambuNDR1ReducedGTFCounts
        File? bambuGTFCounts = bambuTask.bambuGTFCounts
    }
}
