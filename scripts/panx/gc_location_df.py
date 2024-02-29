import pandas as pd
import utils as ut
import pathlib
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_loci_dfs", nargs="+", type=str, required=True)
    parser.add_argument("--in_gc_json", type=str, required=True)
    parser.add_argument("--gid", type=int, required=True)
    parser.add_argument("--out_df", type=str, required=True)
    return parser.parse_args()


def split_genecluster_locus(loc_tag):
    loc_tag = loc_tag.split("_")
    iso = "_".join(loc_tag[:2])
    loc = "_".join(loc_tag[2:])
    return iso, loc


def parse_loci_dfs(loci_dfs_fnames):
    loci_location_dfs = {}
    for loci_df_fname in loci_dfs_fnames:
        iso = pathlib.Path(loci_df_fname).stem
        loci_location_dfs[iso] = pd.read_csv(loci_df_fname, index_col=0)
    return loci_location_dfs


if __name__ == "__main__":

    args = parse_args()

    # load geneCluster.json
    jc = ut.load_genecluster_json(args.in_gc_json)

    # select jc entry for the gene of interest
    gid = args.gid
    jc = jc.loc[gid]

    # load loci position dataframe for each isolate
    loc_dfs = parse_loci_dfs(args.in_loci_dfs)

    # pairs of (isolate, locus) for each gene cluster
    loci = [split_genecluster_locus(loc) for loc in jc["locus"].split()]

    # create a dataframe with the loci positions
    gdf = []
    for iso, loc in loci:
        df = loc_dfs[iso]
        entry = df.loc[loc].copy()
        entry["iso"] = iso
        gdf.append(entry)
    gdf = pd.concat(gdf, axis=1).T
    gdf.drop(columns=["strand", "length"], inplace=True)
    gdf.rename(columns={"start": "beg"}, inplace=True)
    gdf.index.name = "id"
    gdf.to_csv(args.out_df)
