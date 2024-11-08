# mash distance -> csv
# core genome distance (from alignment) -> csv
# block distance (P/A and length of private sequence) -> csv
# merge csv and start analysis from this


rule DST_corealignment:
    input:
        fa=rules.PG_corealignment.output.fa,
        json=rules.PG_corealignment.output.json,
    output:
        "results/{dset}/distances/coredivergence-{opt}.csv",
    params:
        n_cons=lambda wildcards, input: extract_value(input.json, "n. consensus"),
        label="core_div_naive",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/distances/alignment_distance.py \
            --aln {input.fa} --n_cons {params.n_cons} --csv {output} \
            --label "{params.label}"
        """


rule DST_filtered_corealignment:
    input:
        fa=rules.PG_filtered_corealignment.output.fa,
        json=rules.PG_filtered_corealignment.output.info_size,
    output:
        "results/{dset}/distances/coredivergence-filtered-{opt}.csv",
    params:
        n_cons=lambda wildcards, input: extract_value(
            input.json, "polished aln consensus"
        ),
        label="core_div_filtered",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/distances/alignment_distance.py \
            --aln {input.fa} --n_cons {params.n_cons} --csv {output} \
            --label "{params.label}"
        """


rule DST_mash:
    input:
        fa=lambda w: expand(rules.gbk_to_fa.output.fa, acc=dset_chrom_accnums[w.dset]),
    output:
        "results/{dset}/distances/mash_dist.csv",
    conda:
        "../conda_env/tree_inference.yml"
    shell:
        """
        mash triangle {input} > {output}.tmp
        python3 scripts/utils/mash_triangle_to_csv.py \
            --mash_tri {output}.tmp --csv {output}
        rm {output}.tmp
        """


rule DST_pangraph:
    input:
        rules.PG_polish.output,
    output:
        "results/{dset}/distances/pangraph{opt}.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/distances/pangraph_pairwise_distance.py \
            --pangraph {input} --csv {output}
        """


rule DST_pangraph_nodupl:
    input:
        rules.PG_polish.output,
    output:
        "results/{dset}/distances/pangraph{opt}-nodupl.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/distances/pangraph_pairwise_distance.py \
            --pangraph {input} --csv {output} \
            --exclude_dupl
        """


rule DST_edge:
    input:
        rules.PG_polish.output,
    output:
        "results/{dset}/distances/pangraph{opt}_edge_distance.csv",
    params:
        len_thr=100,
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/distances/edge_distance.py \
            --pan {input} --csv {output} --len_thr {params.len_thr}
        """


rule DST_blocks:
    input:
        rules.PG_polish.output,
    output:
        "results/{dset}/distances/pangraph{opt}_block_distance.csv",
    params:
        len_thr=100,
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/distances/block_distance.py \
            --pan {input} --csv {output} --len_thr {params.len_thr}
        """


rule DST_merge:
    input:
        aln=rules.DST_corealignment.output,
        aln_flt=rules.DST_filtered_corealignment.output,
        mash=rules.DST_mash.output,
        pan=rules.DST_pangraph.output,
        edge=rules.DST_edge.output,
        block=rules.DST_blocks.output,
    output:
        "results/{dset}/distances/summary-{opt}.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/distances/merge_dist.py \
            --dfs_in {input.aln} {input.mash} {input.pan} \
                     {input.aln_flt} {input.edge} {input.block} \
            --csv {output}
        """


rule DST_all:
    input:
        expand(rules.DST_merge.output, dset=dset_names, opt=kernel_opts),
        expand(rules.DST_pangraph_nodupl.output, dset=dset_names, opt=kernel_opts),
