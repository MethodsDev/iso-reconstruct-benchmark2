
ref_guided:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/talon TALON-runner.py --output_prefix testSample  --genome /workspace/data/SIRVs1-7.genome.fa --gtf /workspace/data/updated_SIRV_annotation.gtf  --bam /workspace/data/sim.fasta.mm2.bam --output_dir /workspace




clean:
	sudo rm -rf ./testSample* ./talon_tmp ./tmp_label_reads
