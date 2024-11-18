res/%_alignment_ungapped.fa: res/%_alignment_ungapped.fa.gz
	gzip -dc $< > $@

res/%_subtree.nwk: res/%_alignment_restricted.fa
	conda run -n fasttree fasttree -nt -gtr $< > $@ 2> log_$*.txt

res/%_subtree_full.nwk: res/%_alignment_ungapped.fa
	conda run -n fasttree fasttree -nt -gtr $< > $@ 2> log_$*.txt

res/%_subtree_optimized.nwk: res/%_alignment_restricted.fa res/%_alignment_lengths.csv res/%_subtree.nwk
	conda run -n treetime python 03_refine_tree.py \
		--tree_in res/$*_subtree.nwk \
		--aln res/$*_alignment_restricted.fa \
		--lengths res/$*_alignment_lengths.csv \
		--tree_out $@


res/%_subtree_optimized_full.nwk: res/%_alignment_ungapped.fa res/%_alignment_lengths.csv res/%_subtree_full.nwk
	conda run -n treetime python 03_refine_tree.py \
		--tree_in res/$*_subtree_full.nwk \
		--aln res/$*_alignment_ungapped.fa \
		--lengths res/$*_alignment_lengths.csv \
		--tree_out $@

res/maintree.nwk:
	cp ../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk $@

res/%_tanglegram: res/maintree.nwk res/%_subtree_optimized.nwk
	treeknit $^ -o $@ --auspice-view

res/%_tanglegram_full: res/maintree.nwk res/%_subtree_optimized_full.nwk
	treeknit $^ -o $@ --auspice-view

all: res/merger_tanglegram res/flanking_tanglegram

all_full: res/merger_tanglegram_full res/flanking_tanglegram_full