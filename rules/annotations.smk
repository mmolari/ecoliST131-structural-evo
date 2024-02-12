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


rule GM_preformat:
    input:
        rules.GM_summary.output,
    output:
        "results/{dset}/annotations/loc/genomad.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/annotations/genomad_df_preformat.py \
            --input_df {input} \
            --output_df {output}
        """


rule IF_annotate:
    input:
        fa=rules.gbk_to_fa.output.fa,
    output:
        if_dir=directory("data/integron_finder/Results_Integron_Finder_{acc}"),
        I="data/integron_finder/Results_Integron_Finder_{acc}/{acc}.integrons",
        S="data/integron_finder/Results_Integron_Finder_{acc}/{acc}.summary",
    params:
        outdir="data/integron_finder",
    conda:
        "../conda_env/integron_finder.yml"
    shell:
        """
        integron_finder \
            --outdir {params.outdir} \
            --circ \
            --local-max \
            --mute \
            --cpu 1 \
            --pdf \
            --func-annot \
            --promoter-attI \
            {input.fa}
        """


rule IF_summary:
    input:
        S=lambda w: expand(rules.IF_annotate.output.S, acc=dset_chrom_accnums[w.dset]),
        I=lambda w: expand(rules.IF_annotate.output.I, acc=dset_chrom_accnums[w.dset]),
    output:
        i_summ="results/{dset}/annotations/integron_finder/integron_summary.tsv",
        i_ann="results/{dset}/annotations/integron_finder/integron_annotations.tsv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        cat {input.S} | grep -v '^#' | tsv-uniq > {output.i_summ}
        cat {input.I} | grep -v '^#' | tsv-uniq > {output.i_ann}
        """


rule IF_preformat:
    input:
        rules.IF_summary.output.i_ann,
    output:
        "results/{dset}/annotations/loc/integronfinder.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/annotations/integronfinder_df_preformat.py \
            --input_df {input} \
            --output_df {output}
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


rule ISEScan_preformat:
    input:
        rules.ISEScan_summary.output,
    output:
        "results/{dset}/annotations/loc/ISEScan.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/annotations/IS_df_preformat.py \
            --input_df {input} \
            --output_df {output}
        """


rule Dfinder_models:
    output:
        directory("data/defensefinder_models"),
    conda:
        "../conda_env/defensefinder.yml"
    shell:
        """
        defense-finder update --models-dir {output}
        """


rule Dfinder_find:
    input:
        fa=rules.gbk_to_fa.output.fa,
        mod=rules.Dfinder_models.output,
    output:
        a=directory("data/defense_finder/{acc}"),
        g="data/defense_finder/{acc}/{acc}_defense_finder_genes.tsv",
        s="data/defense_finder/{acc}/{acc}_defense_finder_systems.tsv",
    conda:
        "../conda_env/defensefinder.yml"
    shell:
        """
        defense-finder run \
            -o {output.a} \
            --models-dir {input.mod} \
            {input.fa}
        """


rule Dfinder_summary:
    input:
        s=lambda w: expand(rules.Dfinder_find.output.s, acc=dset_chrom_accnums[w.dset]),
        g=lambda w: expand(rules.Dfinder_find.output.g, acc=dset_chrom_accnums[w.dset]),
    output:
        s="results/{dset}/annotations/defense_finder/systems_summary.tsv",
        g="results/{dset}/annotations/defense_finder/genes_summary.tsv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        tsv-append -H {input.s} > {output.s}
        tsv-append -H {input.g} > {output.g}
        """


zero_based_tools = {
    "ISEScan": False,
    "genomad": False,
    "integronfinder": False,  # to be checked
}


rule AN_assign_positions:
    input:
        el="results/{dset}/annotations/loc/{tool}.csv",
        j_pos=rules.BJ_extract_joints_pos.output.pos,
        iso_len=rules.PG_genome_lengths.output,
    output:
        "results/{dset}/annotations/junct_pos_{opt}/{tool}_{K}.csv",
    conda:
        "../conda_env/bioinfo.yml"
    params:
        zero_based=lambda w: "--zero_based" if zero_based_tools[w.tool] else "",
        random=lambda w: "--random" if w.K == "rand" else "",
    shell:
        """
        python3 scripts/annotations/assing_junction.py \
            --iso_len {input.iso_len} \
            --junction_pos_json {input.j_pos} \
            --element_pos_df {input.el} \
            --output_pos {output} \
            {params.zero_based} \
            {params.random}
        """


rule AN_all:
    input:
        expand(rules.GM_summary.output, dset="ST131_ABC"),
        expand(rules.IF_summary.output, dset="ST131_ABC"),
        expand(rules.Dfinder_summary.output, dset="ST131_ABC"),
        expand(
            rules.AN_assign_positions.output,
            dset="ST131_ABC",
            tool=["ISEScan", "genomad", "integronfinder"],
            opt=["asm20-100-5"],
            K=["rand", "real"],
        ),
