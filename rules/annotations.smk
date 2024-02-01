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
        d=directory("data/genomad/{acc}"),
        s="data/genomad/{acc}/{acc}_summary/{acc}_virus_summary.tsv",
    conda:
        "../conda_env/genomad.yml"
    shell:
        """
        genomad end-to-end {input.fa} {output.d} {input.db} \
            --cleanup \
            --threads 4
        """


rule GM_summary:
    input:
        lambda w: expand(rules.GM_run.output.s, acc=dset_chrom_accnums[w.dset]),
    output:
        "results/{dset}/annotations/genomad/prophage_summary.tsv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        tsv-append -H {input} > {output}
        """


rule ISEScan_run:
    input:
        fa=rules.gbk_to_fa.output.fa,
    output:
        d=directory("data/ISEScan/{acc}"),
        s="data/ISEScan/{acc}/fa/{acc}.fa.tsv",
    conda:
        "../conda_env/isescan.yml"
    shell:
        """
        isescan.py --seqfile {input.fa} --output {output.d} --nthread 6
        """


rule ISEScan_summary:
    input:
        lambda w: expand(rules.ISEScan_run.output.s, acc=dset_chrom_accnums[w.dset]),
    output:
        "results/{dset}/annotations/isescan/is_summary.tsv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        tsv-append -H {input} > {output}
        """


rule AN_all:
    input:
        expand(rules.GM_summary.output, dset="ST131_ABC"),
        expand(rules.ISEScan_summary.output, dset="ST131_ABC"),
