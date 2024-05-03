version 1.0

task MandalorionTask {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        String datasetName
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 128
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/mandalorion@sha256:9a2dd74d2a716ed59784b75f64ddf43e451e59e0afb31dfe40176eed4a2460cf"
        File monitoringScript = "gs://ctat_genome_libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

#        samtools bam2fq ~{inputBAM} > samtools.bam2fq.fastq
#        sed -n '1~4s/^@/>/p;2~4p' samtools.bam2fq.fastq > samtools.bam2fq.fasta

        python3 /usr/local/src/Mandalorion/utils/removePolyA_nonDirectionalInput.py -i ~{inputBAM} -o samtools.bam2fq.noPolyA.fasta -t 1,1

        /usr/local/src/Mandalorion/minimap2/misc/paftools.js gff2bed ~{referenceAnnotation} > anno.bed
        /usr/local/src/Mandalorion/minimap2/minimap2 -G 400k --secondary=no -uf -ax splice:hq --cs=long --junc-bed anno.bed -t ~{numThreads} ~{referenceGenome} samtools.bam2fq.noPolyA.fasta > samtools.view.sam
        

        rm samtools.bam2fq.fastq
        rm samtools.bam2fq.fasta
        rm samtools.bam2fq.noPolyA.fasta
        rm anno.bed
        samtools view -bS samtools.view.sam > "~{inputBAM}.bam"
        rm samtools.view.sam
    >>>

    output {
        File bamfile = "~{inputBAM}.bam"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}

workflow Mandalorion {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File referenceAnnotation
        String datasetName
    }

    call MandalorionTask {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation = referenceAnnotation,
            datasetName = datasetName,
    }

    output {
        File bamfile = MandalorionTask.bamfile
        File monitoringLog = MandalorionTask.monitoringLog
    }
}
