rule PL_resistance:
    input:
        tabs=lambda w: expand(
            rules.RG_abricate.output.tab,
            iso=plasmid_accnums[w.dset],
            allow_missing=True,
        ),
    output:
        txt="results/{dset}/plasmids/resistance/{database}_summary.txt",
    conda:
        "../conda_env/abricate.yml"
    shell:
        """
        abricate --summary {input.tabs} > {output.txt}
        """


rule PL_all:
    input:
        expand(
            rules.PL_resistance.output.txt,
            dset=plasmid_accnums.keys(),
            database=["card", "ncbi"],
        ),
