rule PX_download_repo:
    output:
        directory("data/panX"),
    shell:
        "git clone git@github.com:neherlab/pan-genome-analysis.git {output}"


rule PX_link_gbk_file:
    input:
        lambda w: expand(rules.download_gbk.output, acc=dset_chrom_accnums[w.dset]),
    output:
        directory("data/panX/data/{dset}/input_GenBank"),
    shell:
        """
        mkdir -p {output}
        for f in {input}; do
            ln -s $f {output}/$(basename $f)
        done
        """
