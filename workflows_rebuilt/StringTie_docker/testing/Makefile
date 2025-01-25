


ref_free:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/stringtie:latest stringtie-runner.py --genome /workspace/data/SIRVs1-7.genome.fa --output_prefix /workspace/testSample.ref_free --bam /workspace/data/sim.fasta.mm2.bam


ref_guided:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/stringtie:latest stringtie-runner.py --genome /workspace/data/SIRVs1-7.genome.fa --output_prefix /workspace/testSample.ref_guided --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --bam /workspace/data/sim.fasta.mm2.bam 


quant_only:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/stringtie:latest  stringtie-runner.py --genome /workspace/data/SIRVs1-7.genome.fa --output_prefix /workspace/testSample.quant_only --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --bam /workspace/data/sim.fasta.mm2.bam --quant_only 
