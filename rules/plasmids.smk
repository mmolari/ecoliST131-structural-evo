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


rule PL_mob_typing:
    input:
        fa=rules.gbk_to_fa.output.fa,
    output:
        txt="data/mob/{acc}_typing.txt",
    conda:
        "../conda_env/mob_suite.yml"
    shell:
        """
        mob_typer --infile {input.fa} --out_file {output.txt}
        """


rule PL_mob_typing_summary:
    input:
        txts=lambda w: expand(
            rules.PL_mob_typing.output.txt,
            acc=plasmid_accnums[w.dset],
        ),
    output:
        tsv="results/{dset}/plasmids/mob/plasmid_summary.tsv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        header=$(head -n 1 {input.txts[0]})
        echo "$header" > {output.tsv}
        tail -n +2 -q {input.txts} >> {output.tsv}
        """


rule PL_alleles_concat:
    input:
        lambda w: expand(
            rules.QC_alleles_assign.output,
            allele=w.allele,
            acc=sorted(plasmid_accnums[w.dset]),
        ),
    output:
        "results/{dset}/plasmids/alleles/{allele}.tsv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        echo "iso\tlocus\tmatch\tallele\tcoverage\tsim\tmatches\taln_L\tallele_L" > {output}
        cat {input} >> {output}
        """


rule PL_alleles_summary:
    input:
        dfs=expand(
            rules.PL_alleles_concat.output,
            allele=config["plsm_alleles"],
            allow_missing=True,
        ),
    output:
        csv="results/{dset}/plasmids/alleles_summary.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/assembly_qc/alleles_summary.py \
            --in_dfs {input.dfs} \
            --summary_df {output.csv}
        """


rule PL_all:
    input:
        expand(
            rules.PL_join_resistance.output,
            dset=plasmid_accnums.keys(),
            database=["card", "ncbi"],
        ),
        expand(
            rules.PL_mob_typing_summary.output,
            dset=plasmid_accnums.keys(),
        ),
        expand(
            rules.PL_alleles_summary.output,
            dset=plasmid_accnums.keys(),
        ),
