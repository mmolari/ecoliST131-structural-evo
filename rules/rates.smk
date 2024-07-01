rule RT_coldspot_df:
    input:
        jdf=rules.BJ_junct_stats.output.stats,
        gm=expand(
            rules.AN_assign_positions.output,
            tool="genomad",
            K="real",
            allow_missing=True,
        ),
        ig=expand(
            rules.AN_assign_positions.output,
            tool="integronfinder",
            K="real",
            allow_missing=True,
        ),
        IS=expand(
            rules.AN_assign_positions.output,
            tool="ISEScan",
            K="real",
            allow_missing=True,
        ),
        df=expand(
            rules.AN_assign_positions.output,
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
        edf="results/{dset}/rates/{opt}/terminal_branches.csv",
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


rule RT_events_df:
    input:
        idf=rules.RT_internal_coldspots_mugration.output.jdf,
        tdf=rules.RT_terminal_coldspots.output.cdf,
        bdf=rules.RT_internal_coldspots_mugration.output.bdf,
    output:
        df="results/{dset}/rates/{opt}/merged_events_df.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/rates/merge_event_df.py \
            --internal_df {input.idf} \
            --terminal_df {input.tdf} \
            --branch_df {input.bdf} \
            --out_df {output.df}
        """


rule FG_junctions_ann:
    input:
        jdf=rules.BJ_junct_stats.output.stats,
        epg=rules.BJ_extract_pangenome_info.output.info,
        gm=expand(
            rules.AN_assign_positions.output,
            tool="genomad",
            K="real",
            allow_missing=True,
        ),
        ig=expand(
            rules.AN_assign_positions.output,
            tool="integronfinder",
            K="real",
            allow_missing=True,
        ),
        IS=expand(
            rules.AN_assign_positions.output,
            tool="ISEScan",
            K="real",
            allow_missing=True,
        ),
        df=expand(
            rules.AN_assign_positions.output,
            tool="defensefinder",
            K="real",
            allow_missing=True,
        ),
    output:
        ff=directory("figs/{dset}/{opt}/junctions/"),
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/junctions.py \
            --junctions_stats {input.jdf} \
            --edge_pangenome {input.epg} \
            --ann_gm {input.gm} \
            --ann_if {input.ig} \
            --ann_is {input.IS} \
            --ann_df {input.df} \
            --out_fld {output.ff}
        """


rule FG_rates:
    input:
        csdf=rules.RT_coldspot_df.output.df,
        ev=rules.RT_events_df.output.df,
        tjdf=rules.RT_terminal_coldspots.output.cdf,
        ijdf=rules.RT_internal_coldspots_mugration.output.jdf,
        tbdf=rules.RT_terminal_coldspots.output.edf,
        ibdf=rules.RT_internal_coldspots_mugration.output.bdf,
        alni=rules.PG_filtered_corealignment.output.info_size,
    output:
        ff=directory("figs/{dset}/{opt}/rates"),
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/rates.py \
            --coldspots_df {input.csdf} \
            --events_df {input.ev} \
            --terminal_j_df {input.tjdf} \
            --internal_j_df {input.ijdf} \
            --terminal_b_df {input.tbdf} \
            --internal_b_df {input.ibdf} \
            --alignment_info {input.alni} \
            --out_fld {output.ff}
        """


rule RT_all:
    input:
        expand(
            rules.RT_events_df.output,
            dset=dset_names,
            opt=kernel_opts,
        ),
        expand(
            rules.FG_junctions_ann.output,
            dset=dset_names,
            opt=kernel_opts,
        ),
        expand(
            rules.FG_rates.output,
            dset=dset_names,
            opt=kernel_opts,
        ),
