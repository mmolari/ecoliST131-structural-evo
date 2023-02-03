
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


rule PG_all:
    input:
        expand(rules.PG_polish.output, dset=datasets.keys(), opt=kernel_opt.keys()),
        expand(rules.PG_export.output, dset=datasets.keys(), opt=kernel_opt.keys()),