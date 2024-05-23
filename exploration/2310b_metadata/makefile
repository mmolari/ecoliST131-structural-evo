figs1:
	mkdir figs1
	conda run -n pangenome python3 metadata.py \
		--metadata_csv ../../results/ST131/metadata.csv \
		--alleles_csv ../../results/ST131/assembly_qc/alleles_summary.csv \
		--coregenome_tree ../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk \
		--outdir figs1

figs2:
	mkdir figs2
	conda run -n pangenome python3 metadata.py \
		--metadata_csv ../../results/ST131_full/metadata.csv \
		--alleles_csv ../../results/ST131_full/assembly_qc/alleles_summary.csv \
		--coregenome_tree ../../results/ST131_full/pangraph/asm20-100-5-filtered-coretree.nwk \
		--outdir figs2

clean:
	rm -r figs1
	rm -r figs2

all:
	make figs1
	make figs2