rule FG_assembly_qc:
    input:
        qc_csv=rules.QC_summary.output,
    output:
        fig="figs/{dset}/assembly_qc.pdf",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/assembly_qc.py \
            --csv {input.qc_csv} --fig {output.fig}
        """


rule FG_metadata:
    input:
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
        meta=rules.metadata_preprocess.output,
        alls=rules.QC_alleles_summary.output.csv,
    output:
        fld=directory("figs/{dset}/{opt}/metadata"),
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/metadata_plots.py \
            --metadata_csv {input.meta} \
            --alleles_csv {input.alls} \
            --coregenome_tree {input.tree} \
            --outdir {output.fld}
        """


rule FG_resistance:
    input:
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
        ncbi=expand(rules.RG_summary.output.txt, database="ncbi", allow_missing=True),
        card=expand(rules.RG_summary.output.txt, database="card", allow_missing=True),
    output:
        fld=directory("figs/{dset}/{opt}/resistance"),
    params:
        thr=config["resistance"]["id_threshold"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/resistance_plots.py \
            --coregenome_tree {input.tree} \
            --ncbi_df {input.ncbi} \
            --card_df {input.card} \
            --id_threshold {params.thr} \
            --outdir {output.fld}
        """


rule FG_plasmid_resistance:
    input:
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
        chrm=rules.RG_summary.output.txt,
        plsm=rules.PL_join_resistance.output.csv,
    output:
        fig="figs/{dset}/{opt}/plasmids/plasmid_resistance_{database}.pdf",
    params:
        thr=config["resistance"]["id_threshold"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/plasmid_resistance.py \
            --coregenome_tree {input.tree} \
            --chr_df {input.chrm} \
            --pls_df {input.plsm} \
            --id_threshold {params.thr} \
            --fig {output.fig}
        """


rule FG_plasmid_mlst:
    input:
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
        df=rules.PL_alleles_summary.output.csv,
        json=lambda w: dsets_config[w.dset]["plasmids"],
    output:
        fig="figs/{dset}/{opt}/plasmids/plasmid_mlst.pdf",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/plasmid_mlst.py \
            --coregenome_tree {input.tree} \
            --mlst_df {input.df} \
            --plsm_json {input.json} \
            --fig {output.fig}
        """


rule FG_recombination_filter:
    input:
        info_idxs=rules.PG_filtered_corealignment.output.info_idxs,
        info_size=rules.PG_filtered_corealignment.output.info_size,
    output:
        full="figs/{dset}/{opt}/pangraph/corealn_remove_recombination_full.pdf",
        reduced="figs/{dset}/{opt}/pangraph/corealn_remove_recombination_reduced.pdf",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/recombination_filter.py \
            --idxs {input.info_idxs} \
            --size {input.info_size} \
            --fig_full {output.full} \
            --fig_reduced {output.reduced}
        """


rule FG_tree_summary:
    input:
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
        mtd=rules.metadata_preprocess.output,
        alle=rules.QC_alleles_summary.output.csv,
        pl_all=rules.PL_alleles_summary.output.csv,
        plsm=lambda w: dsets_config[w.dset]["plasmids"],
        r_pls=expand(
            rules.PL_join_resistance.output.csv, database="card", allow_missing=True
        ),
        r_chr=expand(rules.RG_summary.output.txt, database="card", allow_missing=True),
    output:
        fig=directory("figs/{dset}/{opt}/tree_summary"),
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/tree.py \
            --tree {input.tree} \
            --metadata {input.mtd} \
            --alleles {input.alle} \
            --plasmids_alleles {input.pl_all} \
            --plasmids_json {input.plsm} \
            --resistance_pls {input.r_pls} \
            --resistance_chr {input.r_chr} \
            --fig_fld {output.fig}
        """


rule FG_homoplasies:
    input:
        tree=rules.PG_coregenome_tree.output.nwk,
        aln=rules.PG_corealignment.output.fa,
        aln_info=rules.PG_corealignment.output.json,
        filt_tree=rules.PG_filtered_coregenome_tree.output.nwk,
        filt_aln=rules.PG_filtered_corealignment.output.fa,
        filt_aln_info=rules.PG_filtered_corealignment.output.info_size,
    output:
        hist_fig="figs/{dset}/{opt}/homoplasies/hist.pdf",
        tree_fig="figs/{dset}/{opt}/homoplasies/tree.pdf",
        homoplasies_fig="figs/{dset}/{opt}/homoplasies/homoplasies.pdf",
    conda:
        "../conda_env/tree_inference.yml"
    shell:
        """
        python3 scripts/figs/homoplasies.py \
            --tree {input.tree} \
            --aln {input.aln} \
            --aln_info {input.aln_info} \
            --filt_tree {input.filt_tree} \
            --filt_aln {input.filt_aln} \
            --filt_aln_info {input.filt_aln_info} \
            --hist_fig {output.hist_fig} \
            --tree_fig {output.tree_fig} \
            --homoplasies_fig {output.homoplasies_fig}
        """


rule FG_block_distr_fig:
    input:
        rules.PG_polish.output,
    output:
        "figs/{dset}/{opt}/pangraph/block_distr.pdf",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/pangraph/plot_block_distr.py \
            --pangraph {input} --fig {output}
        """


rule FG_distances:
    input:
        csv=rules.DST_merge.output,
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
    output:
        directory("figs/{dset}/{opt}/distances"),
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/distances.py \
            --dist_df {input.csv} \
            --tree {input.tree} \
            --fig_fld {output}
        """


rule FG_coresynt:
    input:
        pg=rules.PG_polish.output,
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
    output:
        fg=directory("figs/{dset}/{opt}/coresynt"),
        mg="results/{dset}/pangraph/coresynt-{opt}/mergers.csv",
        bc="results/{dset}/pangraph/coresynt-{opt}/blocks.csv",
    params:
        len_thr=config["backbone-joints"]["len-thr"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/fig_core_synteny.py \
            --pangraph {input.pg} \
            --tree {input.tree} \
            --len_thr {params.len_thr} \
            --mergers {output.mg} \
            --block_colors {output.bc} \
            --fig_fld {output.fg}
        """


rule FG_circle_synteny:
    input:
        pg=rules.PG_polish.output,
        fa=lambda w: expand(
            rules.gbk_to_fa.output, acc=config["datasets"][w.dset]["guide-strain"]
        ),
        bc=rules.FG_coresynt.output.bc,
        mg=rules.FG_coresynt.output.mg,
    output:
        fg="figs/{dset}/{opt}/circle_synteny.png",
    params:
        len_thr=config["backbone-joints"]["len-thr"],
        ref=lambda w: config["datasets"][w.dset]["guide-strain"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/synteny_circle.py \
            --pan {input.pg} \
            --ref {params.ref} \
            --fa {input.fa} \
            --len_thr {params.len_thr} \
            --block_colors {input.bc} \
            --mergers {input.mg} \
            --fig {output.fg}
        """


rule FG_junctions_survey:
    input:
        df=rules.BJ_extract_joints_df.output.dfl,
    output:
        directory("figs/{dset}/{opt}/junctions_survey"),
    params:
        len_thr=config["backbone-joints"]["len-thr"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/fig_junctions_survey.py \
            --joints_df {input.df} \
            --fig_fld {output}
        """


rule FG_junctions_stats:
    input:
        jdf=rules.BJ_junct_stats.output.stats,
        epg=rules.BJ_extract_pangenome_info.output.info,
    output:
        ff=directory("figs/{dset}/{opt}/junctions_stats/"),
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/junctions_simple.py \
            --junctions_stats {input.jdf} \
            --edge_pangenome {input.epg} \
            --out_fld {output.ff}
        """


rule FG_all:
    input:
        expand(rules.FG_assembly_qc.output, dset=dset_names),
        expand(rules.FG_metadata.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_homoplasies.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_recombination_filter.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_resistance.output, dset=dset_names, opt=kernel_opts),
        expand(
            rules.FG_plasmid_resistance.output,
            dset=plasmid_dsets.keys(),
            opt=kernel_opts,
            database=["ncbi", "card"],
        ),
        expand(rules.FG_plasmid_mlst.output, dset=plasmid_dsets.keys(), opt=kernel_opts),
        expand(rules.FG_tree_summary.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_block_distr_fig.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_distances.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_coresynt.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_circle_synteny.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_junctions_survey.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_junctions_stats.output, dset=dset_names, opt=kernel_opts),
