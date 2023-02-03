
kernel_opt = config["pangraph"]["kernel-options"]
datasets = config["datasets"]

wildcard_constraints:
    opt=f"({'|'.join(kernel_opt.keys())})",
    dset=f"({'|'.join(datasets.keys())})",


rule PG_build:
    input:
        fa=lambda w: expand(rules.gbk_to_fa.output.fa, acc=datasets[w.dset]),
    output:
        "results/pangraph/{dset}/{opt}.json",
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
        "results/pangraph/{dset}/{opt}-polished.json",
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
        directory("results/pangraph/{dset}/{opt}_export"),
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
        expand(rules.PG_polish.output, dset=datasets.keys(), opt=kernel.keys()),
        expand(rules.PG_export.output, dset=datasets.keys(), opt=kernel.keys()),