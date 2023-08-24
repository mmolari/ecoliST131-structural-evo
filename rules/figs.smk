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


rule FG_all:
    input:
        expand(rules.FG_assembly_qc.output, dset=datasets.keys()),
        expand(rules.FG_metadata.output, dset=datasets.keys(), opt=kernel_opt.keys()),
