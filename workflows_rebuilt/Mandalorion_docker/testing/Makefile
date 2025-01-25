

ref_guided:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/mandalorion Mandalorian-runner.py --output_prefix /workspace/ref_guided --genome /workspace/data/SIRVs1-7.genome.fa --bam /workspace/data/sim.fasta.mm2.bam --gtf /workspace/data/updated_SIRV_annotation.gtf --delayTime 1


ref_free:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/mandalorion Mandalorian-runner.py --output_prefix /workspace/ref_free --genome /workspace/data/SIRVs1-7.genome.fa --bam /workspace/data/sim.fasta.mm2.bam --delayTime 1


