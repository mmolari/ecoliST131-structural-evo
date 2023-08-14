rule BJ_find_edges:
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


rule BJ_extract_joints_pos:
    input:
        pan=rules.PG_polish.output,
        edges=rules.BJ_find_edges.output.edges,
    output:
        pos="results/{dset}/backbone_joints/{opt}/joints_pos.json",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/extract_joints_positions.py \
            --pangraph {input.pan} \
            --edges {input.edges} \
            --positions {output.pos}
        """


rule B_extract_raw_paths:
    input:
        pan=rules.PG_polish.output,
        edges=rules.BJ_find_edges.output.edges,
    output:
        paths="results/{dset}/backbone_joints/{opt}/raw_paths.json",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/extract_raw_paths.py \
            --pangraph {input.pan} \
            --edges {input.edges} \
            --paths_json {output.paths}
        """


rule BJ_all:
    input:
        expand(
            rules.BJ_extract_joints_pos.output,
            dset=datasets.keys(),
            opt=kernel_opt.keys(),
        ),
        expand(
            rules.B_extract_raw_paths.output,
            dset=datasets.keys(),
            opt=kernel_opt.keys(),
        ),
