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


rule FG_recombination_filter:
    input:
        info_idxs=rules.PG_filtered_corealignment.output.info_idxs,
        info_size=rules.PG_filtered_corealignment.output.info_size,
    output:
        full="figs/{dset}/pangraph/corealn_remove_recombination/{opt}_full.pdf",
        reduced="figs/{dset}/pangraph/corealn_remove_recombination/{opt}_reduced.pdf",
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
        "figs/{dset}/pangraph/{opt}_block_distr.pdf",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/pangraph/plot_block_distr.py \
            --pangraph {input} --fig {output}
        """


rule FG_all:
    input:
        expand(rules.FG_assembly_qc.output, dset=dset_names),
        expand(rules.FG_metadata.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_homoplasies.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_recombination_filter.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_resistance.output, dset=dset_names, opt=kernel_opts),
        expand(rules.FG_block_distr_fig.output, dset=dset_names, opt=kernel_opts),
