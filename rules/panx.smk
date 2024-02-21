rule PX_download_repo:
    output:
        "data/panX/README.md",
    shell:
        "git clone https://github.com/neherlab/pan-genome-analysis.git $(dirname {output})"


rule PX_link_gbk_file:
    input:
        gbk=lambda w: expand(rules.download_gbk.output, acc=dset_chrom_accnums[w.dset]),
        rp=rules.PX_download_repo.output,
    output:
        directory("data/panX/data/{dset}/input_GenBank"),
    shell:
        """
        mkdir -p {output}
        for f in {input.gbk}; do
            ln -s $(realpath $f) {output}/$(basename $f)
        done
        """


rule PX_run:
    input:
        rp=rules.PX_download_repo.output,
        d=rules.PX_link_gbk_file.output,
    output:
        gcj="data/panX/data/{dset}/vis/geneCluster.json",
        tree="data/panX/data/{dset}/vis/strain_tree.nwk",
        gcd=directory("data/panX/data/{dset}/vis/geneCluster"),
        log="data/panX/data/{dset}/log.txt",
        err="data/panX/data/{dset}/err.txt",
    conda:
        "../conda_env/panX-environment.yml"
    shell:
        """
        python $(dirname {input.rp})/panX.py \
            -fn $(dirname {input.d}) \
            -t 10 \
            -sl {wildcards.dset} \
            > {output.log} 2> {output.err}
        """


rule PX_loci_df:
    input:
        gbk=rules.download_gbk.output,
    output:
        "data/loci_df/{acc}.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        "python scripts/panx/loci_df.py --in_gbk {input.gbk} --out_df {output}"


rule PX_all:
    input:
        expand(rules.PX_run.output, dset="ST131_ABC"),
        expand(rules.PX_loci_df.output, acc=dset_chrom_accnums["ST131_ABC"]),
