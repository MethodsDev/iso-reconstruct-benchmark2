
test: ref_guided quant_only

wdl_test: wdl_ref_guided wdl_quant_only

ref_guided:
	docker run --rm -i -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/flair:latest flair-runner.py --genome /workspace/data/SIRVs1-7.genome.fa --bam /workspace/data/sim.fasta.mm2.bam --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --output_prefix /workspace/testSample.ref_guided


quant_only:
	docker run --rm -i -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/flair:latest flair-runner.py --genome /workspace/data/SIRVs1-7.genome.fa --bam /workspace/data/sim.fasta.mm2.bam --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --output_prefix /workspace/testSample.quant_only --quant_only


wdl_ref_guided:
	miniwdl run ../../Flair.wdl -i wdl_inputs/ref_guided.inputs.json


wdl_quant_only:
	miniwdl run ../../Flair.wdl -i wdl_inputs/quant_only.inputs.json
