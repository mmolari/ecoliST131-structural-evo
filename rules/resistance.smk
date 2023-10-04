rule RG_abricate:
    input:
        gbk="data/gbk/{iso}.gbk",
    output:
        tab="data/resistance/{database}/{iso}.tab",
    conda:
        "../conda_env/abricate.yml"
    shell:
        """
        abricate --db {wildcards.database} {input.gbk} > {output.tab}
        """


rule RG_summary:
    input:
        tabs=lambda w: expand(
            rules.RG_abricate.output.tab,
            iso=dset_chrom_accnums[w.dset],
            allow_missing=True,
        ),
    output:
        txt="results/{dset}/resistance/{database}_summary.txt",
    conda:
        "../conda_env/abricate.yml"
    shell:
        """
        abricate --summary {input.tabs} > {output.txt}
        """


rule RG_all:
    input:
        expand(rules.RG_summary.output.txt, dset=dset_names, database=["card", "ncbi"]),
