test: quant_only ref_guided ref_free


quant_only:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant IsoQuant-runner.py  --genome /workspace/data/SIRVs1-7.genome.fa --bam /workspace/data/sim.fasta.mm2.bam --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --data_type pacbio_ccs --output_prefix /workspace/quant_only --quant_only


ref_guided:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant IsoQuant-runner.py  --genome /workspace/data/SIRVs1-7.genome.fa --bam /workspace/data/sim.fasta.mm2.bam --gtf /workspace/data/SIRVs1-7.annot.reduced.gtf --data_type pacbio_ccs --output_prefix /workspace/ref_guided


ref_free:
	docker run --rm -it -v `pwd`:/workspace us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/isoquant IsoQuant-runner.py  --genome /workspace/data/SIRVs1-7.genome.fa --bam /workspace/data/sim.fasta.mm2.bam  --data_type pacbio_ccs --output_prefix /workspace/ref_free



clean:
	sudo rm -rf ./quant_only.* ./ref_free.* ./ref_guided.*
