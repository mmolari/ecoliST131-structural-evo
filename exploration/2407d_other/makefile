joint_graph_fld := results/ST131_ABC/backbone_joints/asm20-100-5/joints_pangraph
tree := results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk
jpos := results/ST131_ABC/backbone_joints/asm20-100-5/joints_pos.json
gbks := data/gbk
ann := results/ST131_ABC/annotations/junct_pos_asm20-100-5
gen_lengths := results/ST131_ABC/pangraph/genome_lengths.csv
hh_genes := exploration/2405a_hochhauser/data/hochhauser.csv


data:
	mkdir -p data

data/hh_genes.csv: data
	ln -sf ../../../$(hh_genes) data/hh_genes.csv

data/tree.nwk: data
	ln -sf ../../../$(tree) data/tree.nwk

data/joints_pos.json: data
	ln -sf ../../../$(jpos) data/joints_pos.json

data/gbk: data
	ln -sf ../../../$(gbks) data/gbk

data/genome_lengths.csv: data
	ln -sf ../../../$(gen_lengths) data/genome_lengths.csv

data/ann: data
	mkdir -p data/ann

data/ann/%.csv: data/ann
	ln -sf ../../../../$(ann)/$*_real.csv data/ann/$*.csv

res/%/joint_graph.json:
	mkdir -p res/$*
	ln -sf ../../../../$(joint_graph_fld)/$*.json $@

res/%/gbk_annotations.csv res/%/tool_annotations.csv: data/genome_lengths.csv data/gbk data/ann/defensefinder.csv data/ann/genomad.csv data/ann/integronfinder.csv data/ann/ISEScan.csv data/joints_pos.json res/%/joint_graph.json
	python3 01_extract_annotations.py --hs $*

res/%/gene_1.png res/%/gene_2.png: res/%/gbk_annotations.csv res/%/tool_annotations.csv res/%/joint_graph.json data/tree.nwk data/hh_genes.csv
	python3 02_draw_all_genes.py --hs $*

LAB := \
	AQOVQHWFSE_f__KIVDDRSTJR_f \
	CUSJYRGUQV_r__XDZOAXPBLI_f \
	CRALSVTFTC_f__VKGJGCHOXW_f \
	FHAPZVUFCF_r__HTGOZRIQGS_r \
	FZZJSRRGUR_r__WNLYKTZGKW_r \
	GAKJONQISG_r__UJBGATQZSS_r \
	KUDTBWWBJE_f__YSEMYTMKWD_f \
	HLDRILYIZH_r__YSKQIEWBCG_r \
	HLZGFRDRSW_r__SYQZDBITDR_r \
	MRDJEEINGI_f__PHTUBOWZWA_r
TGTS := $(foreach lab,$(LAB),res/$(lab)/gene_1.png)

.PHONY: all
all:
	$(MAKE) $(TGTS)
