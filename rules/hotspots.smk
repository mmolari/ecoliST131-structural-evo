import itertools as itt


checkpoint HP_select_hotspots:
    input:
        stats=rules.BJ_junct_stats.output.stats,
        pang=rules.BJ_extract_pangenome_info.output.info,
    output:
        hs="results/{dset}/hotspots/{opt}/hotspots.csv",
        fig="results/{dset}/hotspots/{opt}/hotspots.png",
    params:
        min_len=config["hotspots"]["min_length"],
        min_npaths=config["hotspots"]["min_npaths"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/hotspots/select_spots.py \
            --j_stats {input.stats} \
            --j_pangenome {input.pang} \
            --out_csv {output.hs} \
            --out_fig {output.fig} \
            --min_len {params.min_len} \
            --min_paths {params.min_npaths}
        """


def read_hotspots(csv_fname):
    idxs = []
    with open(csv_fname, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            row = line.split(",")
            if len(row) > 1:
                idxs.append(row[0].strip())
    return idxs


rule HP_hotspot_stats:
    input:
        hs_pan=rules.BJ_pangraph.output.pan,
        dst=rules.DST_merge.output,
    output:
        csv="results/{dset}/hotspots/{opt}/hotspot_stats/{edge}.csv",
        fig="results/{dset}/hotspots/{opt}/hotspot_stats/{edge}.png",
    params:
        mbl=config["hotspots"]["min_block_length"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/hotspots/extract_stats.py \
            --pan {input.hs_pan} \
            --pair_dist {input.dst} \
            --min_block_len {params.mbl} \
            --out_csv {output.csv} \
            --out_fig {output.fig}
        """


# def all_hotspots_stats(wildcards):
#     hotspot_file = checkpoints.HP_select_hotspots.get(**wildcards).output["hs"]
#     hotspots = read_hotspots(hotspot_file)
#     return expand(rules.HP_hotspot_stats.output.csv, edge=hotspots, **wildcards)


def all_hotspots_stats(wildcards):
    files = []
    for dset, opt in itt.product(dset_names, kernel_opts):
        wc = {"dset": dset, "opt": opt}
        hotspot_file = checkpoints.HP_select_hotspots.get(**wc).output["hs"]
        hotspots = read_hotspots(hotspot_file)
        files += expand(rules.HP_hotspot_stats.output.csv, edge=hotspots, **wc)
    return files


rule HP_HH_locate:
    input:
        hh=config["hotspots"]["hochhauser_SI"],
        gbk=lambda w: expand(rules.download_gbk.output, acc=dset_chrom_accnums[w.dset]),
    output:
        csv="results/{dset}/hotspots/hochhauser_coregenes_position.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/hotspots/hochhauser_locate.py \
            --hochhauser {input.hh} \
            --out {output.csv} \
            --gbk {input.gbk}
        """


rule HP_all:
    input:
        all_hotspots_stats,
        expand(rules.HP_HH_locate.output, dset=dset_names),
