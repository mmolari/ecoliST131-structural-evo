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


rule AN_all:
    input:
        expand(rules.GM_summary.output, dset="ST131_ABC"),
