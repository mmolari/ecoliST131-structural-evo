rule download_gbk:
    output:
        "data/gbk/{acc}.gbk",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        ncbi-acc-download {wildcards.acc} -e all -F genbank
        mv {wildcards.acc}.gbk {output}
        """


rule gbk_to_fa:
    input:
        gbk=rules.download_gbk.output,
    output:
        fa="data/fa/{acc}.fa",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/utils/gbk_to_fa.py --gbk {input.gbk} --fa {output.fa}
        """


rule metadata_preprocess:
    input:
        lambda w: dsets_config[w.dset]["metadata"],
    output:
        "results/{dset}/metadata.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/utils/metadata_cleanup.py \
            --metadata_in {input} \
            --metadata_out {output}
        """
