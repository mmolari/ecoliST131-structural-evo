data/aln.fa:
	ln -s ../../../results/ST131_ABC/pangraph/asm20-100-5-alignment/full_corealignment.fa data/aln.fa

data/tree.nwk:
	ln -s ../../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk data/tree.nwk

data/dates.csv:
	python3 00_dates.py

treetime: data/aln.fa data/tree.nwk data/dates.csv
	treetime --aln data/aln.fa --tree data/tree.nwk --dates data/dates.csv --outdir results

clean:
	rm -rf results
	rm -f data/aln.fa data/tree.nwk data/dates.csv