# load accession numbers of excluded isolates
excluded = {k: [] for k in datasets.keys()}
for k, fname in config["excluded"].items():
    with open(fname, "r") as f:
        acc_nums = f.readlines()
    acc_nums = [an.strip() for an in acc_nums]
    acc_nums = [an for an in acc_nums if len(an) > 0]
    excluded[k] = acc_nums


rule QC_busco_download:
    output:
        directory("data/busco_models"),
    conda:
        "../conda_env/busco.yml"
    shell:
        """
        busco --download enterobacterales_odb10
        mv busco_downloads {output}
        """


# busco --download prokaryota --download_path {output}


rule QC_busco_run:
    input:
        fa=rules.gbk_to_fa.output.fa,
        mod=rules.QC_busco_download.output,
    output:
        directory("results/{dset}/assembly_qc/busco/{acc}"),
    conda:
        "../conda_env/busco.yml"
    shell:
        """
        busco -i {input.fa} \
            --offline \
            -l enterobacterales_odb10 \
            --download_path {input.mod} \
            -m genome \
            -o {output}
        """


# --auto-lineage-prok \


rule QC_summary:
    input:
        lambda w: expand(
            rules.QC_busco_run.output,
            acc=sorted(datasets[w.dset] + excluded[w.dset]),
            allow_missing=True,
        ),
    output:
        "results/{dset}/assembly_qc/summary.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/assembly_qc/busco_summary.py \
            --fld {input} \
            --csv {output}
        """


rule QC_all:
    input:
        expand(rules.QC_summary.output, dset=datasets.keys()),
