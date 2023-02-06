# mash distance -> csv
# core genome distance (from alignment) -> csv
# block distance (P/A and length of private sequence) -> csv
# merge csv and start analysis from this


rule DST_corealignment:
    input:
        fa=rules.PG_reduced_corealignment.output.fa,
        json=rules.PG_reduced_corealignment.output.json,
    output:
        "results/{dset}/distances/{opt}-coredivergence.csv",
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
        "results/{dset}/distances/mash.csv",
    conda:
        "../conda_env/fasttree.yml"
    shell:
        """
        mash triangle {input} > {output}.tmp
        python3 scripts/utils/mash_triangle_to_csv.py \
            --mash_tri {output}.tmp --csv {output}
        rm {output}.tmp
        """


rule DST_all:
    input:
        expand(
            rules.DST_corealignment.output, dset=datasets.keys(), opt=kernel_opt.keys()
        ),
        expand(rules.DST_mash.output, dset=datasets.keys()),
