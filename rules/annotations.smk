rule GM_download_db:
    output:
        db=directory("data/genomad_db"),
    conda:
        "../conda_env/genomad.yml"
    shell:
        """
        genomad download-database data/
        """


rule GM_run:
    input:
        db=rules.GM_download_db.output.db,
        fa=rules.gbk_to_fa.output.fa,
    output:
        "data/genomad/{acc}",
    conda:
        "../conda_env/genomad.yml"
    shell:
        """
        genomad end-to-end {input.fa} {output} {input.db} \
            --cleanup \
            --threads 4 \
            --splits 0
        """


rule AN_all:
    input:
        [
            expand(rules.GM_run.output, acc=dset_chrom_accnums[k])
            for k in dset_chrom_accnums.keys()
        ],
