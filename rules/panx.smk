rule PX_download_repo:
    output:
        directory("data/panX"),
    shell:
        "git clone git@github.com:neherlab/pan-genome-analysis.git {output}"


rule PX_link_gbk_file:
    input:
        gbk=lambda w: expand(rules.download_gbk.output, acc=dset_chrom_accnums[w.dset]),
        rp=rules.PX_download_repo.output,
    output:
        main=directory("data/panX/data/{dset}"),
        data=directory("data/panX/data/{dset}/input_GenBank"),
    shell:
        """
        mkdir -p {output.data}
        for f in {input.gbk}; do
            ln -s $f {output.data}/$(basename $f)
        done
        """

rule PX_run:
    input:
        d=rules.PX_link_gbk_file.output.main,
    output:
        gcj="data/panX/data/{dset}/vis/geneCluster.json",
        tree="data/panX/data/{dset}/vis/strain_tree.nwk",
        gcd=directory("data/panX/data/{dset}/vis/geneCluster"),
        log="data/panX/data/{dset}/log.txt",
        err="data/panX/data/{dset}/err.txt",
    conda:
        "../conda_env/panX-environment.yml",
    shell:
        """
        python {input.rp}/panX.py \
            -fn {input.d} \
            -t 10 \
            -sl {wildcards.dset} \
            > {output.log} 2> {output.err}
        """
```