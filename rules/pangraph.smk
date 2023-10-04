rule PG_build:
    input:
        fa=lambda w: expand(rules.gbk_to_fa.output.fa, acc=dset_chrom_accnums[w.dset]),
    output:
        "results/{dset}/pangraph/{opt}.json",
    params:
        opt=lambda w: config["pangraph"]["kernel-options"][w.opt],
    shell:
        """
        export JULIA_NUM_THREADS=8
        pangraph build {params.opt} {input.fa} > {output}
        """


rule PG_polish:
    input:
        rules.PG_build.output,
    output:
        "results/{dset}/pangraph/{opt}-polished.json",
    params:
        opt=config["pangraph"]["polish-options"],
    conda:
        "../conda_env/pangraph.yml"
    shell:
        """
        export JULIA_NUM_THREADS=8
        pangraph polish {params.opt} {input} > {output}
        """


rule PG_export:
    input:
        rules.PG_polish.output,
    output:
        directory("results/{dset}/pangraph/{opt}_export"),
    conda:
        "../conda_env/pangraph.yml"
    shell:
        """
        pangraph export \
            --no-duplications \
            --output-directory {output} \
            --prefix export \
            {input}
        """


rule PG_block_distr_fig:
    input:
        rules.PG_polish.output,
    output:
        "figs/{dset}/pangraph/{opt}_block_distr.pdf",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/pangraph/plot_block_distr.py \
            --pangraph {input} --fig {output}
        """


rule PG_corealignment:
    input:
        rules.PG_polish.output,
    output:
        fa="results/{dset}/pangraph/{opt}-alignment/corealignment.fa",
        json="results/{dset}/pangraph/{opt}-alignment/corealignment_info.json",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/pangraph/reduced_core_alignment.py \
            --pangraph {input} --fasta_aln {output.fa} --info {output.json}
        """


rule PG_filtered_corealignment:
    input:
        rules.PG_polish.output,
    output:
        fa="results/{dset}/pangraph/{opt}-alignment/filtered_corealignment.fa",
        info_size="results/{dset}/pangraph/{opt}-alignment/filtered_corealignment_info_size.json",
        info_idxs="results/{dset}/pangraph/{opt}-alignment/filtered_corealignment_info_idxs.json",
    conda:
        "../conda_env/bioinfo.yml"
    params:
        window=1000,
        max_nsnps=3,
        guide_strain=lambda w: dsets_config[w.dset]["guide-strain"],
    shell:
        """
        python3 scripts/pangraph/corealn_without_recombination.py \
            --pangraph {input} \
            --guide_strain {params.guide_strain} \
            --window {params.window} \
            --max_nsnps {params.max_nsnps} \
            --fasta_aln {output.fa} \
            --info_size {output.info_size} \
            --info_idxs {output.info_idxs}
        """


def extract_value(json_file, key):
    with open(json_file, "r") as f:
        data = json.load(f)
    return data[key]


rule PG_coregenome_tree:
    input:
        fa=rules.PG_corealignment.output.fa,
        json=rules.PG_corealignment.output.json,
    output:
        nwk="results/{dset}/pangraph/{opt}-coretree.nwk",
    conda:
        "../conda_env/tree_inference.yml"
    params:
        n_cons=lambda wildcards, input: extract_value(input.json, "n. consensus"),
    shell:
        """
        fasttree -gtr -nt {input.fa} > {output.nwk}.tmp
        python3 scripts/pangraph/refine_coretree_treetime.py \
            --tree_in {output.nwk}.tmp \
            --n_consensus {params.n_cons} \
            --aln {input.fa} \
            --tree_out {output.nwk}
        rm {output.nwk}.tmp
        """


rule PG_filtered_coregenome_tree:
    input:
        fa=rules.PG_filtered_corealignment.output.fa,
        json=rules.PG_filtered_corealignment.output.info_size,
    output:
        nwk="results/{dset}/pangraph/{opt}-filtered-coretree.nwk",
    conda:
        "../conda_env/tree_inference.yml"
    params:
        n_cons=lambda wildcards, input: extract_value(
            input.json, "polished aln consensus"
        ),
    shell:
        """
        fasttree -gtr -nt {input.fa} > {output.nwk}.tmp
        python3 scripts/pangraph/refine_coretree_treetime.py \
            --tree_in {output.nwk}.tmp \
            --n_consensus {params.n_cons} \
            --aln {input.fa} \
            --tree_out {output.nwk}
        rm {output.nwk}.tmp
        """


rule PG_all:
    input:
        expand(rules.PG_block_distr_fig.output, dset=dset_names, opt=kernel_opts),
        expand(rules.PG_coregenome_tree.output, dset=dset_names, opt=kernel_opts),
        expand(
            rules.PG_filtered_coregenome_tree.output, dset=dset_names, opt=kernel_opts
        ),
        expand(rules.PG_export.output, dset=dset_names, opt=kernel_opts),
