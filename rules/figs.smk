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
    output:
        hist="figs/{dset}/{opt}/metadata/hist.pdf",
        tree="figs/{dset}/{opt}/metadata/tree.pdf",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/metadata_tree.py \
            --tree {input.tree} \
            --metadata {input.meta} \
            --hist_fig {output.hist} \
            --tree_fig {output.tree}
        """


rule FG_homoplasies:
    input:
        tree=rules.PG_coregenome_tree.output.nwk,
        aln=rules.PG_corealignment.output.fa,
        aln_info=rules.PG_corealignment.output.json,
        filt_tree=rules.PG_filtered_coregenome_tree.output.nwk,
        filt_aln=rules.PG_filtered_corealignment.output.fa,
        filt_aln_info=rules.PG_filtered_corealignment.output.json,
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


rule FG_all:
    input:
        expand(rules.FG_assembly_qc.output, dset=datasets.keys()),
        expand(rules.FG_metadata.output, dset=datasets.keys(), opt=kernel_opt.keys()),
        expand(rules.FG_homoplasies.output, dset=datasets.keys(), opt=kernel_opt.keys()),
        expand(
            rules.PG_filtered_corealignment.output.plot_reduced,
            dset=datasets.keys(),
            opt=kernel_opt.keys(),
        ),
