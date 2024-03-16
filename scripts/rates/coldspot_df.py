import pandas as pd
import numpy as np
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--junctions_stats", type=str)
    parser.add_argument("--ann_genomad", type=str)
    parser.add_argument("--ann_integrons", type=str)
    parser.add_argument("--ann_isescan", type=str)
    parser.add_argument("--ann_defensefinder", type=str)
    parser.add_argument("--out_df", type=str)

    return parser.parse_args()


def load_mge_annotations(args):
    return {
        "genomad": pd.read_csv(args.ann_genomad, index_col=0),
        "integrons": pd.read_csv(args.ann_integrons, index_col=0),
        "isescan": pd.read_csv(args.ann_isescan, index_col=0),
        "defensefinder": pd.read_csv(args.ann_defensefinder, index_col=0),
    }


def add_mge_annotations(df, mge_annotations):
    for k, df2 in mge_annotations.items():
        df2 = df2["junction"].value_counts()
        df[f"{k}"] = df.index.map(df2)
        df[f"{k}"] = df[f"{k}"].fillna(0)
        df[f"{k}"] = df[f"{k}"].astype(int)

    return df


if __name__ == "__main__":
    args = parse_args()

    jdf = pd.read_csv(args.junctions_stats, index_col=0)

    # add annotations
    mges = load_mge_annotations(args)
    jdf = add_mge_annotations(jdf, mges)

    # keep only coldspots
    mask = jdf["n_categories"] == 2
    # remove non-backbone
    mask &= jdf["n_iso"] == jdf["n_iso"].max()
    jdf = jdf[mask]

    assert not np.any(jdf["transitive"]), "transitive junctions should be absent"

    # drop irrelevant fields
    jdf = jdf.drop(
        columns=[
            "n_iso",
            "has_dupl",
            "n_categories",
            "cat_entropy",
            "n_nodes",
            "min_length",
            "max_length",
            "mean_length",
            "core_left_length",
            "core_right_length",
            "transitive",
            "nonempty_acc_len",
        ]
    )

    # sort values
    jdf = jdf.sort_values("majority_category", ascending=False)

    # save to file
    jdf.to_csv(args.out_df)
