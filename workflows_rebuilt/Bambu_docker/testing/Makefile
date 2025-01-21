quant_only:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/bambu:latest bambu-runner.Rscript --bam /workspace/data/sim.fasta.mm2.bam --genome /workspace/data/SIRVs1-7.genome.fa --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf  --output_prefix /workspace/testSample --quant_only


ref_free:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/bambu:latest bambu-runner.Rscript --bam /workspace/data/sim.fasta.mm2.bam --genome /workspace/data/SIRVs1-7.genome.fa  --output_prefix /workspace/testSample


ref_guided:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/bambu:latest  bambu-runner.Rscript --bam /workspace/data/sim.fasta.mm2.bam --genome /workspace/data/SIRVs1-7.genome.fa --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf  --output_prefix /workspace/testSample


