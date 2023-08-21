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
        directory("results/{dset}/assembly_qc/{acc}"),
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
        lambda w:expand(rules.QC_busco_run.output, acc=datasets[w.dset]),
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
        [expand(rules.QC_busco_run.output, dset=dset, acc=acc) for dset, acc in datasets.items()],
