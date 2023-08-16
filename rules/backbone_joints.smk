import itertools as itt

BJ_config = config["backbone-joints"]


checkpoint BJ_find_edges:
    input:
        pan=rules.PG_polish.output,
    output:
        edges="results/{dset}/backbone_joints/{opt}/core_edges.csv",
    params:
        len_thr=BJ_config["len-thr"],
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


rule BJ_extract_raw_paths:
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


rule BJ_extract_joint_sequence:
    input:
        genomes=lambda w: expand(rules.gbk_to_fa.output.fa, acc=datasets[w.dset]),
        pos=rules.BJ_extract_joints_pos.output.pos,
    output:
        seq="results/{dset}/backbone_joints/{opt}/joints_seqs/{edge}.fa",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/extract_joint_sequence.py \
            --edge_pos {input.pos} \
            --edge {wildcards.edge} \
            --out_fa {output.seq} \
            --genomes {input.genomes}
        """


rule BJ_pangraph:
    input:
        seq=rules.BJ_extract_joint_sequence.output.seq,
    output:
        pan="results/{dset}/backbone_joints/{opt}/joints_pangraph/{edge}.json",
    params:
        opt_build=BJ_config["build-opt"],
        opt_polish=BJ_config["polish-opt"],
    conda:
        "../conda_env/pangraph.yml"
    shell:
        """
        export JULIA_NUM_THREADS=1
        pangraph build {params.opt_build} {input.seq} | \
        pangraph polish {params.opt_polish} > {output.pan}
        """


rule BJ_mash_dist:
    input:
        seq=rules.BJ_extract_joint_sequence.output.seq,
    output:
        dist="results/{dset}/backbone_joints/{opt}/dist/mash/{edge}.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        mash triangle {input.seq} > {output.dist}.tmp
        python3 scripts/utils/mash_triangle_to_csv.py \
            --mash_tri {output.dist}.tmp --csv {output.dist}
        rm {output.dist}.tmp
        """


rule BJ_plot_linear_repr:
    input:
        pan=rules.BJ_pangraph.output.pan,
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
    output:
        fig="figs/{dset}/backbone_joints/{opt}/joints_linear_plot/{edge}.png",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/plot_junction_categories.py \
            --pangraph {input.pan} \
            --tree {input.tree} \
            --fig {output.fig}
        """


def BJ_all_joints_outputs(wildcards):
    files = []
    for dset, opt in itt.product(datasets.keys(), kernel_opt.keys()):
        wc = {"dset": dset, "opt": opt}

        # define list of edges
        edge_file = checkpoints.BJ_find_edges.get(**wc).output["edges"]
        with open(edge_file, "r") as f:
            edges = []
            for line in f.readlines():
                edge, n = line.strip().split(",")
                if n == str(len(datasets[dset])):
                    edges.append(edge)

        # add desired output files
        files += expand(rules.BJ_plot_linear_repr.output.fig, edge=edges, **wc)
        files += expand(rules.BJ_mash_dist.output.dist, edge=edges, **wc)
    return files


rule BJ_all:
    input:
        expand(
            rules.BJ_extract_joints_pos.output,
            dset=datasets.keys(),
            opt=kernel_opt.keys(),
        ),
        BJ_all_joints_outputs,
