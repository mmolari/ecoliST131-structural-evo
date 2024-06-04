data/aln.fa:
	ln -s ../../../results/ST131_ABC/pangraph/asm20-100-5-alignment/filtered_corealignment.fa data/aln.fa

data/tree.nwk:
	ln -s ../../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk data/tree.nwk

data/dates.csv:
	python3 00_dates.py

treetime_1: data/aln.fa data/tree.nwk data/dates.csv
	treetime --aln data/aln.fa \
		--tree data/tree.nwk \
		--dates data/dates.csv \
		--outdir results_1 \
		--sequence-length 2427416

treetime_2: data/aln.fa data/tree.nwk data/dates.csv
	treetime --aln data/aln.fa \
		--tree data/tree.nwk \
		--dates data/dates.csv \
		--outdir results_2 \
		--clock-rate 2.46e-7 \
		--sequence-length 2427416

treetime_3: data/aln.fa data/tree.nwk data/dates.csv
	treetime --aln data/aln.fa \
		--tree data/tree.nwk \
		--dates data/dates.csv \
		--outdir results_3 \
		--clock-rate 4.39e-7 \
		--sequence-length 2427416

all: treetime_1 treetime_2 treetime_3

clean:
	rm -rf results_1
	rm -rf results_2
	rm -rf results_3
	rm -f data/aln.fa data/tree.nwk data/dates.csv