rule RT_coldspot_df:
    input:
        jdf=rules.BJ_junct_stats.output.stats,
        gm=expand(
            "results/{dset}/annotations/junct_pos_{opt}/{tool}_{K}.csv",
            tool="genomad",
            K="real",
            allow_missing=True,
        ),
        ig=expand(
            "results/{dset}/annotations/junct_pos_{opt}/{tool}_{K}.csv",
            tool="integronfinder",
            K="real",
            allow_missing=True,
        ),
        IS=expand(
            "results/{dset}/annotations/junct_pos_{opt}/{tool}_{K}.csv",
            tool="ISEScan",
            K="real",
            allow_missing=True,
        ),
        df=expand(
            "results/{dset}/annotations/junct_pos_{opt}/{tool}_{K}.csv",
            tool="defensefinder",
            K="real",
            allow_missing=True,
        ),
    output:
        df="results/{dset}/rates/{opt}/coldspot_df.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/rates/coldspot_df.py \
            --junctions_stats {input.jdf} \
            --ann_genomad {input.gm} \
            --ann_integrons {input.ig} \
            --ann_isescan {input.IS} \
            --ann_defensefinder {input.df} \
            --out_df {output.df}
        """


rule RT_all:
    input:
        expand(
            rules.RT_coldspot_df.output,
            dset=config["datasets"],
            opt=config["pangraph"]["kernel-options"],
        ),
