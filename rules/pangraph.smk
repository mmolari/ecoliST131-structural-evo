kernel_opt = config["pangraph"]["kernel-options"]
datasets_fnames = config["datasets"]

# read accession numbers from dataset files
datasets = {}
for k, fname in datasets_fnames.items():
    with open(fname, "r") as f:
        acc_nums = f.readlines()
    acc_nums = [an.strip() for an in acc_nums]
    acc_nums = [an for an in acc_nums if len(an) > 0]
    datasets[k] = acc_nums


wildcard_constraints:
    opt=f"({'|'.join(kernel_opt.keys())})",
    dset=f"({'|'.join(datasets.keys())})",


rule PG_build:
    input:
        fa=lambda w: expand(rules.gbk_to_fa.output.fa, acc=datasets[w.dset]),
    output:
        "results/{dset}/pangraph/{opt}.json",
    params:
        opt=lambda w: kernel_opt[w.opt],
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


rule PG_reduced_corealignment:
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


rule PG_all:
    input:
        expand(
            rules.PG_block_distr_fig.output,
            dset=datasets.keys(),
            opt=kernel_opt.keys(),
        ),
        expand(
            rules.PG_reduced_corealignment.output,
            dset=datasets.keys(),
            opt=kernel_opt.keys(),
        ),
        expand(rules.PG_export.output, dset=datasets.keys(), opt=kernel_opt.keys()),
