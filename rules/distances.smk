# mash distance -> csv
# core genome distance (from alignment) -> csv
# block distance (P/A and length of private sequence) -> csv
# merge csv and start analysis from this


rule DST_corealignment:
    input:
        fa=rules.PG_reduced_corealignment.output.fa,
        json=rules.PG_reduced_corealignment.output.json,
    output:
        "results/{dset}/distances/coredivergence-{opt}.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/distances/alignment_distance.py \
            --aln {input.fa} --info {input.json} --csv {output}
        """


rule DST_mash:
    input:
        fa=lambda w: expand(rules.gbk_to_fa.output.fa, acc=datasets[w.dset]),
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


rule DST_merge:
    input:
        aln=rules.DST_corealignment.output,
        mash=rules.DST_mash.output,
        pan=rules.DST_pangraph.output,
    output:
        "results/{dset}/distances/summary-{opt}.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/distances/merge_dist.py \
            --dfs_in {input.aln} {input.mash} {input.pan} \
            --csv {output}
        """


rule DST_all:
    input:
        expand(rules.DST_merge.output, dset=datasets.keys(), opt=kernel_opt.keys()),
