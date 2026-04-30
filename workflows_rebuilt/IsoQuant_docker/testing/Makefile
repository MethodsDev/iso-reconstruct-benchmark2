test: quant_only ref_guided ref_free

wdl_test: wdl_quant_only wdl_ref_guided wdl_ref_free

quant_only:
	docker run --rm -i -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:latest IsoQuant-runner.py  --genome /workspace/data/SIRVs1-7.genome.fa --bam /workspace/data/sim.fasta.mm2.bam --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --data_type pacbio_ccs --output_prefix /workspace/quant_only --quant_only


ref_guided:
	docker run --rm -i -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:latest IsoQuant-runner.py  --genome /workspace/data/SIRVs1-7.genome.fa --bam /workspace/data/sim.fasta.mm2.bam --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --data_type pacbio_ccs --output_prefix /workspace/ref_guided


ref_free:
	docker run --rm -i -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant:latest IsoQuant-runner.py  --genome /workspace/data/SIRVs1-7.genome.fa --bam /workspace/data/sim.fasta.mm2.bam  --data_type pacbio_ccs --output_prefix /workspace/ref_free

wdl_quant_only:
	miniwdl run ../../IsoQuant.wdl -i wdl_inputs/quant_only.inputs.json


wdl_ref_guided:
	miniwdl run ../../IsoQuant.wdl -i wdl_inputs/ref_guided.inputs.json


wdl_ref_free:
	miniwdl run ../../IsoQuant.wdl -i wdl_inputs/ref_free.inputs.json



clean:
	sudo rm -rf ./quant_only.* ./ref_free.* ./ref_guided.* ./202* ./call-*
