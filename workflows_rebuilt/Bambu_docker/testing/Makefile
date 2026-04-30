test: quant_only ref_free ref_guided

wdl_test: wdl_quant_only wdl_ref_guided wdl_ref_free

quant_only:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/bambu:latest bambu-runner.Rscript --bam /workspace/data/sim.fasta.mm2.bam --genome /workspace/data/SIRVs1-7.genome.fa --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf  --output_prefix /workspace/testSample --quant_only


ref_free:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/bambu:latest bambu-runner.Rscript --bam /workspace/data/sim.fasta.mm2.bam --genome /workspace/data/SIRVs1-7.genome.fa  --output_prefix /workspace/testSample


ref_guided:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/bambu:latest  bambu-runner.Rscript --bam /workspace/data/sim.fasta.mm2.bam --genome /workspace/data/SIRVs1-7.genome.fa --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf  --output_prefix /workspace/testSample


wdl_quant_only:
	miniwdl run ../../Bambu.wdl -i wdl_inputs/quant_only.inputs.json


wdl_ref_guided:
	miniwdl run ../../Bambu.wdl -i wdl_inputs/ref_guided.inputs.json


wdl_ref_free:
	miniwdl run ../../Bambu.wdl -i wdl_inputs/ref_free.inputs.json


clean:
	sudo rm -rf ./202* ./call-* ./_LAST
