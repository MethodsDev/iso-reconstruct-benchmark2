test: ref_guided quant_only

wdl_test: wdl_ref_guided wdl_quant_only

ref_guided:
	docker run --rm -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/flames-v2:latest \
	--output_prefix /workspace/testSample.ref_guided \
	--genome /workspace/data/SIRVs1-7.genome.fa \
	--bam /workspace/data/sim.fasta.mm2.bam \
	--gtf /workspace/data/SIRVs1-7.annot.reduced.gtf \
	--data_type pacbio_ccs


quant_only:
	docker run --rm -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/flames-v2:latest \
	--output_prefix /workspace/testSample.quant_only \
	--genome /workspace/data/SIRVs1-7.genome.fa \
	--bam /workspace/data/sim.fasta.mm2.bam \
	--gtf /workspace/data/SIRVs1-7.annot.reduced.gtf \
	--data_type pacbio_ccs \
	--quant_only


wdl_ref_guided:
	miniwdl run ../../FLAMESv2.wdl -i wdl_inputs/ref_guided.inputs.json


wdl_quant_only:
	miniwdl run ../../FLAMESv2.wdl -i wdl_inputs/quant_only.inputs.json


clean:
	sudo rm -rf ./202* ./call-* ./_LAST ./flamesv2_out ./testSample.ref_guided.* ./testSample.quant_only.*
