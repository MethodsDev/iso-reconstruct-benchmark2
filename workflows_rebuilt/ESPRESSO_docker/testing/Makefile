test: ref_guided

wdl_test: wdl_ref_guided

ref_guided:
	docker run --rm -i -v `pwd`:/workspace -v `pwd`/..:/espresso_local us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/espresso:latest python /espresso_local/ESPRESSO-runner.py --output_prefix /workspace/testSample.ref_guided --genome /workspace/data/SIRVs1-7.genome.fa --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --bam /workspace/data/sim.fasta.mm2.bam


wdl_ref_guided:
	miniwdl run ../../ESPRESSO.wdl -i wdl_inputs/ref_guided.inputs.json

clean:
	sudo rm -rf ./202* ./call-* ./_LAST ./testSample.ref_guided.*
