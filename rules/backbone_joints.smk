import itertools as itt


BJ_config = config["backbone-joints"]


def read_edge_count(csv_fname):
    idxs = []
    with open(csv_fname, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            edge, count = line.split(",")
            if int(count) > 1:
                idxs.append(edge.strip())
    return idxs


checkpoint BJ_extract_joints_df:
    input:
        pan=rules.PG_polish.output,
    output:
        dfl="results/{dset}/backbone_joints/{opt}/edge_len.csv",
        dfc="results/{dset}/backbone_joints/{opt}/edge_count.csv",
    params:
        len_thr=BJ_config["len-thr"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/junctions_df_extract.py \
            --pangraph {input.pan} \
            --len_thr {params.len_thr} \
            --df_len {output.dfl} \
            --df_count {output.dfc}
        """


rule BJ_extract_joints_pos:
    input:
        pan=rules.PG_polish.output,
        edges_len=rules.BJ_extract_joints_df.output.dfl,
    output:
        pos="results/{dset}/backbone_joints/{opt}/joints_pos.json",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/extract_joints_positions.py \
            --pangraph {input.pan} \
            --edge_len_df {input.edges_len} \
            --positions {output.pos}
        """


# rule BJ_extract_raw_paths:
#     input:
#         pan=rules.PG_polish.output,
#         edges=rules.BJ_find_edges.output.edges,
#     output:
#         paths="results/{dset}/backbone_joints/{opt}/raw_paths.json",
#     conda:
#         "../conda_env/bioinfo.yml"
#     shell:
#         """
#         python3 scripts/backbone_joints/extract_raw_paths.py \
#             --pangraph {input.pan} \
#             --edges {input.edges} \
#             --paths_json {output.paths}
#         """


rule BJ_extract_pangenome_info:
    input:
        pan=rules.PG_polish.output,
        dfl=rules.BJ_extract_joints_df.output.dfl,
    output:
        info="results/{dset}/backbone_joints/{opt}/edge_pangenome.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/extract_joints_pangenome.py \
            --pangraph {input.pan} \
            --edge_df {input.dfl} \
            --out_edge_len_df {output.info}
        """


rule BJ_extract_joint_sequence:
    input:
        genomes=lambda w: expand(
            rules.gbk_to_fa.output.fa, acc=dset_chrom_accnums[w.dset]
        ),
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
        export JULIA_NUM_THREADS=4
        pangraph build {params.opt_build} {input.seq} | \
        pangraph polish {params.opt_polish} > {output.pan}
        """


def all_junction_pangraphs(wildcards):
    edge_count_file = checkpoints.BJ_extract_joints_df.get(**wildcards).output["dfc"]
    edges = read_edge_count(edge_count_file)
    return expand(rules.BJ_pangraph.output.pan, edge=edges, **wildcards)


rule BJ_junct_stats:
    input:
        pans=all_junction_pangraphs,
    output:
        stats="results/{dset}/backbone_joints/{opt}/junctions_stats.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/junctions_stats.py \
            --junct_pangraphs {input.pans} \
            --df_csv {output.stats}
        """


rule BJ_aln_position:
    input:
        pan=rules.PG_polish.output,
        gen_len=rules.PG_genome_lengths.output,
    output:
        SNPs="results/{dset}/backbone_joints/{opt}/j_aln_pos/SNPs.csv",
        junct_pos="results/{dset}/backbone_joints/{opt}/j_aln_pos/junct_pos.csv",
        aln_pos="results/{dset}/backbone_joints/{opt}/j_aln_pos/aln_pos.csv",
    params:
        guide_strain=lambda w: config["datasets"][w.dset]["guide-strain"],
        thr_len=BJ_config["len-thr"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/junction_aln_positions.py \
            --guide_strain {params.guide_strain} \
            --thr_len {params.thr_len} \
            --genome_lengths {input.gen_len} \
            --pangraph {input.pan} \
            --out_SNPs {output.SNPs} \
            --out_junct_pos {output.junct_pos} \
            --out_aln_pos {output.aln_pos}
        """


# rule BJ_mash_dist:
#     input:
#         seq=rules.BJ_extract_joint_sequence.output.seq,
#     output:
#         dist="results/{dset}/backbone_joints/{opt}/dist/mash/{edge}.csv",
#     conda:
#         "../conda_env/bioinfo.yml"
#     shell:
#         """
#         mash triangle {input.seq} > {output.dist}.tmp
#         python3 scripts/utils/mash_triangle_to_csv.py \
#             --mash_tri {output.dist}.tmp --csv {output.dist}
#         rm {output.dist}.tmp
#         """


# rule BJ_plot_linear_repr:
#     input:
#         pan=rules.BJ_pangraph.output.pan,
#         tree=rules.PG_filtered_coregenome_tree.output.nwk,
#     output:
#         fig="figs/{dset}/backbone_joints/{opt}/joints_linear_plot/{edge}.png",
#     conda:
#         "../conda_env/bioinfo.yml"
#     shell:
#         """
#         python3 scripts/backbone_joints/plot_junction_categories.py \
#             --pangraph {input.pan} \
#             --tree {input.tree} \
#             --fig {output.fig}
#         """


def BJ_all_joints_outputs(wildcards):
    files = []
    for dset, opt in itt.product(dset_names, kernel_opts):
        wc = {"dset": dset, "opt": opt}

        # define list of edges
        edge_count_file = checkpoints.BJ_extract_joints_df.get(**wc).output["dfc"]
        edges = read_edge_count(edge_count_file)

        # add desired output files
        # files += expand(rules.BJ_plot_linear_repr.output.fig, edge=edges, **wc)
        # files += expand(rules.BJ_mash_dist.output.dist, edge=edges, **wc)
        files += expand(rules.BJ_pangraph.output, edge=edges, **wc)

    return files


rule BJ_all:
    input:
        expand(rules.BJ_extract_pangenome_info.output, dset=dset_names, opt=kernel_opts),
        expand(rules.BJ_junct_stats.output, dset=dset_names, opt=kernel_opts),
        expand(rules.BJ_aln_position.output, dset=dset_names, opt=kernel_opts),
        # BJ_all_joints_outputs,
