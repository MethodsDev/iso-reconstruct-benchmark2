test: ref_free ref_guided


ref_free:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoseq IsoSeq-runner.py  --output_prefix /workspace/testSample.ref_free --genome /workspace/data/SIRVs1-7.genome.fa --bam /workspace/data/sim.fasta.mm2.bam 

ref_guided:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoseq IsoSeq-runner.py  --output_prefix /workspace/testSample.ref_guided --genome /workspace/data/SIRVs1-7.genome.fa --bam /workspace/data/sim.fasta.mm2.bam --gtf /workspace/data/updated_SIRV_annotation.gtf


clean:
	sudo rm -rf ./testSample*
