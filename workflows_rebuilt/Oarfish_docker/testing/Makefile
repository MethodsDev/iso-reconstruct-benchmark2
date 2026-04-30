
test: byAlignment byReads

wdl_test: wdl_byAlignment wdl_byReads


byAlignment:
	docker run --rm -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/oarfish Oarfish_runner.py --output_prefix /workspace/testSample --genome /workspace/data/SIRVs1-7.genome.fa --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --fastq /workspace/data/reads.fastq.gz --mode byAlignment


byReads:
	docker run --rm -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/oarfish Oarfish_runner.py --output_prefix /workspace/testSample --genome /workspace/data/SIRVs1-7.genome.fa --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --fastq /workspace/data/reads.fastq.gz --mode byReads


wdl_byAlignment:
	miniwdl run ../../Oarfish.wdl -i wdl_inputs/byAlignment.inputs.json


wdl_byReads:
	miniwdl run ../../Oarfish.wdl -i wdl_inputs/byReads.inputs.json


clean:
	sudo rm -rf ./202* ./call-* ./_LAST
