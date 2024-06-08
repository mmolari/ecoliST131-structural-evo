hs=CIRMBUYJFK_f__CWCCKOQCWZ_r
data_dir=data/$(hs)
res_dir=res/$(hs)
fig_dir=figs/$(hs)
joint_graph=results/ST131_ABC/backbone_joints/asm20-100-5/joints_pangraph/$(hs).json
tree=results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk
jpos=results/ST131_ABC/backbone_joints/asm20-100-5/joints_pos.json
gbks=data/gbk
ann=results/ST131_ABC/annotations/junct_pos_asm20-100-5
gen_lengths=results/ST131_ABC/pangraph/genome_lengths.csv
hh_genes=exploration/2405a_hochhauser/data/hochhauser.csv

shared_infiles:
	mkdir -p data
	ln -s ../../../$(tree) data/tree.nwk
	ln -s ../../../$(jpos) data/joints_pos.json
	ln -s ../../../$(gbks) data/gbk
	ln -s ../../../$(gen_lengths) data/genome_lengths.csv
	ln -s ../../../$(hh_genes) data/hh_genes.csv
	mkdir -p data/ann
	ln -s ../../../../$(ann)/defensefinder_real.csv data/ann/defensefinder.csv
	ln -s ../../../../$(ann)/genomad_real.csv data/ann/genomad.csv
	ln -s ../../../../$(ann)/integronfinder_real.csv data/ann/integronfinder.csv
	ln -s ../../../../$(ann)/ISEScan_real.csv data/ann/ISEScan.csv

infiles:
	mkdir -p $(data_dir)
	mkdir -p $(res_dir)
	mkdir -p $(fig_dir)
	ln -s ../../../../$(joint_graph) $(data_dir)/joint_graph.json
	
path_cats:
	python3 00_path_categories.py --hs $(hs)

extract_ann:
	python3 01_extract_annotations.py --hs $(hs)

pangraph_gfa:
	pangraph export -ell 0 --no-duplications \
		--output-directory $(res_dir)/export \
		$(data_dir)/joint_graph.json

draw_all:
	python3 02_draw_all_paths.py --hs $(hs)

draw_genes_all:
	python3 02_draw_all_genes.py --hs $(hs)

all: infiles path_cats extract_ann pangraph_gfa draw_all draw_genes_all

clean_shared:
	rm -r data

clean:
	rm -rf $(data_dir)
	rm -rf $(res_dir)
	rm -rf $(fig_dir)
