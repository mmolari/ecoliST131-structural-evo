rule GB_input_tsv_file:
    input:
        fas=lambda w: expand(rules.gbk_to_fa.output.fa, acc=dset_chrom_accnums[w.dset]),
    output:
        "results/{dset}/gubbins/input_genomes.tsv",
    shell:
        """
        for f in {input.fas}; do
            echo "$(basename $f .fa)\t$f" >> {output}
        done
        """


rule GB_build_ska_alignment:
    input:
        fas=lambda w: expand(rules.gbk_to_fa.output.fa, acc=dset_chrom_accnums[w.dset]),
        tsv=rules.GB_input_tsv_file.output,
        ref=lambda w: expand(
            rules.gbk_to_fa.output.fa, acc=config["datasets"][w.dset]["guide-strain"]
        ),
    output:
        "results/{dset}/gubbins/ska_aln_{dset}.fa",
    conda:
        "../conda_env/gubbins.yml"
    shell:
        """
        generate_ska_alignment.py \
            --reference {input.ref} \
            --input {input.tsv} \
            --out {output}
        """


rule GB_run_gubbins:
    input:
        aln=rules.GB_build_ska_alignment.output,
    output:
        directory("results/{dset}/gubbins/results"),
    params:
        threads=1,
        tree_builder="fasttree",
        prefix=lambda w: f"gubbins_{w.dset}_",
    conda:
        "../conda_env/gubbins.yml"
    shell:
        """
        run_gubbins.py \
            --threads {params.threads} \
            --tree-builder {params.tree_builder} \
            --prefix {params.prefix} \
            --verbose {input.aln}
        mkdir -p {output}
        mv {params.prefix}* {output}
        """


rule GB_all:
    input:
        expand(rules.GB_run_gubbins.output, dset=dset_names),
