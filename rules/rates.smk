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


rule RT_terminal_coldspots:
    input:
        df=rules.RT_coldspot_df.output.df,
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
        jpg=all_junction_pangraphs,
    output:
        cdf="results/{dset}/rates/{opt}/terminal_coldspot.csv",
        edf="results/{dset}/rates/{opt}/terminal_edges.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/rates/terminal_coldspot_analysis.py \
            --coldspots_df {input.df} \
            --tree {input.tree} \
            --junction_pangraphs {input.jpg} \
            --out_coldspot_df {output.cdf} \
            --out_events_df {output.edf} \
        """


rule RT_all:
    input:
        expand(
            rules.RT_terminal_coldspots.output,
            dset=config["datasets"],
            opt=config["pangraph"]["kernel-options"],
        ),
