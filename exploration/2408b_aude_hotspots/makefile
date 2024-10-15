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

res/%/annotations.pdf: res/%/gbk_annotations.csv res/%/tool_annotations.csv res/%/joint_graph.json data/tree.nwk data/hh_genes.csv
	python3 02_draw_all_genes.py --hs $*

LAB := \
	JJRRWBDVGH_f__TFKJQKVUKX_f \
	KAIPBNCIHR_r__WXCHSHHCDT_f \
	HEDBBNUKLU_r__VSCDQDGRHO_r \
	RKAOKULCFF_f__VFXLTFSKTV_r \
	CIRMBUYJFK_f__CWCCKOQCWZ_r \
	GPKQYOCEJI_r__NKVSUZGURN_f \
	CXKKDGMPSE_f__KBLPANZOCZ_f \
	EJPOGALASQ_f__KUIFCLFQSI_r \
	VZTFXIZVXB_f__YOCIMVGHSL_f \
	XFUKGTZLAV_f__YBWQKVQGZE_f \
	ATPWUNKKID_f__KKPYPKGMXA_f \
	KYQOKYBCOW_r__XXIWNZXZTK_r \
	IFRPFFEGON_r__TCWDRAKLPS_r \
	IHKFSQQUKE_r__KPBYGJHRZJ_f \
	XXVMWZCEKI_r__YUOECYBHUS_r \
	JVNRLCFAVD_f__PLTCZQCVRD_r \
	RYYAQMEJGY_r__ZTHKZYHPIX_f \
	

TGTS := $(foreach lab,$(LAB),res/$(lab)/annotations.pdf)

.SECONDARY:

.PHONY: all
all:
	$(MAKE) $(TGTS)