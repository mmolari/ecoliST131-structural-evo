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


checkpoint PX_run:
    input:
        rp=rules.PX_download_repo.output,
        d=rules.PX_link_gbk_file.output,
    output:
        gcj="data/panX/data/{dset}/vis/geneCluster.json",
        tree="data/panX/data/{dset}/vis/strain_tree.nwk",
        acf="data/panX/data/{dset}/allclusters_final.tsv",
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


rule PX_gc_loc_df:
    message:
        "Processing gene cluster {wildcards.gid}"
    input:
        ldf=lambda w: expand(rules.PX_loci_df.output, acc=dset_chrom_accnums[w.dset]),
        gcj=lambda w: expand(rules.PX_run.output.gcj, dset=w.dset),
    output:
        "results/{dset}/panx/gc_loci/genecl_{gid}.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python scripts/panx/gc_location_df.py \
            --in_loci_dfs {input.ldf} \
            --in_gc_json {input.gcj} \
            --gid {wildcards.gid} \
            --out_df {output}
        """


def get_geneclusters_ids(wildcards):
    with checkpoints.PX_run.get(**wildcards).output["gcj"].open() as f:
        gcj = json.load(f)
    gene_ids = [x["geneId"] for x in gcj]
    return expand(rules.PX_gc_loc_df.output, gid=gene_ids, **wildcards)


rule PX_all:
    input:
        expand(rules.PX_run.output, dset="ST131_ABC"),
        lambda w: get_geneclusters_ids({"dset": "ST131_ABC"}),
