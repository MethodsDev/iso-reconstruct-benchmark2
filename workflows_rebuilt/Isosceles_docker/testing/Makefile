ref_guided:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoceles2:latest --bam /workspace/data/sim.fasta.mm2.bam --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --genome /workspace/data/SIRVs1-7.genome.fa --output_prefix /workspace/testing.refguided

quant_only:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoceles2:latest --bam /workspace/data/sim.fasta.mm2.bam --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --genome /workspace/data/SIRVs1-7.genome.fa --output_prefix /workspace/testing.quant_only --quant_only
