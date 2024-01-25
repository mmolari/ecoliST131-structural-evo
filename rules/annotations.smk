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
        directory("data/genomad/{acc}"),
    conda:
        "../conda_env/genomad.yml"
    shell:
        """
        genomad end-to-end {input.fa} {output} {input.db} \
            --cleanup \
            --threads 4
        """


rule GM_summary:
    input:
        lambda w: expand(rules.GM_run.output, acc=dset_chrom_accnums[w.dset]),
    output:
        "results/{dset}/annotations/genomad/prophage_summary.tsv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        ACC=$(basename {input[0]})
        FNAME="{input[0]}/${{ACC}}_summary/${{ACC}}_virus_summary.tsv"
        head -n 1 $FNAME > {output}

        for input_file in {input}; do
            ACC=$(basename $input_file)
            FNAME="$input_file/${{ACC}}_summary/${{ACC}}_virus_summary.tsv"
            tail -n +2 $FNAME >> {output}
        done
        """


rule IF_annotate:
    input:
        fa=rules.gbk_to_fa.output.fa,
    output:
        if_dir=directory("data/integron_finder/Results_Integron_Finder_{acc}"),
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
            {input.fa}
        """


rule IF_summary:
    input:
        lambda w: expand(
            rules.IF_annotate.output.if_dir, acc=dset_chrom_accnums[w.dset]
        ),
    output:
        i_summ="results/{dset}/annotations/integron_finder/integron_summary.tsv",
        i_ann="results/{dset}/annotations/integron_finder/integron_annotations.tsv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        ACC=$(basename {input[0]})
        ACC=${{ACC#Results_Integron_Finder_}}
        FNAME="{input[0]}/${{ACC}}.summary"
        sed -n '2p' $FNAME > {output.i_summ}

        for input_file in {input}; do
            ACC=$(basename $input_file)
            ACC=${{ACC#Results_Integron_Finder_}}
            FNAME="$input_file/${{ACC}}.summary"
            tail -n +3 $FNAME >> {output.i_summ}
        done

        ACC=$(basename {input[0]})
        ACC=${{ACC#Results_Integron_Finder_}}
        FNAME="{input[0]}/${{ACC}}.integrons"
        sed -n '2p' $FNAME > {output.i_ann}

        for input_file in {input}; do
            ACC=$(basename $input_file)
            ACC=${{ACC#Results_Integron_Finder_}}
            FNAME="$input_file/${{ACC}}.integrons"
            tail -n +3 $FNAME >> {output.i_ann}
        done
        """


rule AN_all:
    input:
        expand(rules.GM_summary.output, dset="ST131_ABC"),
        expand(rules.IF_summary.output, dset="ST131_ABC"),
