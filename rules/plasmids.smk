rule PL_resistance:
    input:
        tabs=lambda w: expand(
            rules.RG_abricate.output.tab,
            iso=plasmid_accnums[w.dset],
            allow_missing=True,
        ),
    output:
        txt="results/{dset}/plasmids/resistance/{database}_summary_plasmid.txt",
    conda:
        "../conda_env/abricate.yml"
    shell:
        """
        abricate --summary {input.tabs} > {output.txt}
        """


rule PL_join_resistance:
    input:
        res=rules.PL_resistance.output.txt,
        pls=lambda w: dsets_config[w.dset]["plasmids"],
    output:
        csv="results/{dset}/plasmids/resistance/{database}_summary_chromosome.csv",
    params:
        thr=config["resistance"]["id_threshold"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/plasmids/join_resistance_df.py \
            --plasmid_json {input.pls} \
            --resistance_tsv {input.res} \
            --threshold_id {params.thr} \
            --out_csv {output.csv}
        """


rule PL_all:
    input:
        expand(
            rules.PL_join_resistance.output.csv,
            dset=plasmid_accnums.keys(),
            database=["card", "ncbi"],
        ),
