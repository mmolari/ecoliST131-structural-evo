res/%_subtree.nwk: res/%_alignment_restricted.fa
	conda run -n fasttree fasttree -nt -gtr $< > $@ 2> log_$*.txt

res/%_subtree_optimized.nwk: res/%_alignment_restricted.fa res/%_alignment_lengths.csv res/%_subtree.nwk
	conda run -n treetime python 03_refine_tree.py \
		--tree_in res/$*_subtree.nwk \
		--aln res/$*_alignment_restricted.fa \
		--lengths res/$*_alignment_lengths.csv \
		--tree_out $@

res/maintree.nwk:
	cp ../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk $@

res/%_tanglegram: res/maintree.nwk res/%_subtree_optimized.nwk
	treeknit $^ -o $@ --auspice-view

all: res/merger_tanglegram res/flanking_tanglegram