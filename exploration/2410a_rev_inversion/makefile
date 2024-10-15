res/subtree.nwk: res/merger_alignment_no_gaps.fa
	conda run -n fasttree fasttree -nt -gtr $< > $@ 2> log.txt

res/subtree_optimized.nwk: res/merger_alignment_no_gaps.fa res/subtree.nwk
	conda run -n treetime python 03_refine_tree.py \
		--tree_in res/subtree.nwk \
		--aln res/merger_alignment_no_gaps.fa \
		--tree_out $@

res/maintree.nwk:
	cp ../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk $@

res/tanglegram: res/maintree.nwk res/subtree_optimized.nwk
	treeknit $^ -o $@ --auspice-view
