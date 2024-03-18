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


rule RT_internal_coldspots_prepare:
    input:
        df=rules.RT_coldspot_df.output.df,
        jpg=all_junction_pangraphs,
    output:
        ABs="results/{dset}/rates/{opt}/AB_states.csv",
        ev="results/{dset}/rates/{opt}/AB_nestedness.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/rates/internal_prepare.py \
            --joints_dset {input.df} \
            --joint_pangraphs {input.jpg} \
            --out_AB_states {output.ABs} \
            --out_event_info {output.ev}
        """


rule RT_internal_coldspots_mugration:
    input:
        ABs=rules.RT_internal_coldspots_prepare.output.ABs,
        ev=rules.RT_internal_coldspots_prepare.output.ev,
        jdf=rules.RT_coldspot_df.output.df,
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
    output:
        jdf="results/{dset}/rates/{opt}/nonsingleton_junct_df.csv",
        bdf="results/{dset}/rates/{opt}/nonsingleton_branch_df.csv",
    conda:
        "../conda_env/treetime.yml"
    shell:
        """
        python3 scripts/rates/internal_mugration.py \
            --AB_states {input.ABs} \
            --AB_nestedness {input.ev} \
            --joints_df {input.jdf} \
            --tree {input.tree} \
            --out_nonsingleton_junct_df {output.jdf} \
            --out_nonsingleton_branch_df {output.bdf}
        """


rule RT_all:
    input:
        expand(
            rules.RT_terminal_coldspots.output,
            dset=config["datasets"],
            opt=config["pangraph"]["kernel-options"],
        ),
        expand(
            rules.RT_internal_coldspots_mugration.output,
            dset=config["datasets"],
            opt=config["pangraph"]["kernel-options"],
        ),
