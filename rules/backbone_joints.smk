rule BJ_block_distr_fig:
    input:
        pan=rules.PG_polish.output,
    output:
        edges="results/{dset}/backbone_joints/{opt}/core_edges.csv",
    params:
        len_thr=config["backbone-joints"]["len-thr"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/find_backbone_edges.py \
            --pangraph {input.pan} \
            --csv {output.edges} \
            --len_thr {params.len_thr}
        """


rule BJ_all:
    input:
        expand(
            rules.BJ_block_distr_fig.output.edges,
            dset=datasets.keys(),
            opt=kernel_opt.keys(),
        ),
