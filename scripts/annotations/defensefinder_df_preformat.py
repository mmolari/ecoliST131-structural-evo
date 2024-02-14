import pandas as pd
import numpy as np
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_genes", required=True)
    parser.add_argument("--input_systems", required=True)
    parser.add_argument("--output_df", required=True)
    return parser.parse_args()


def system_beg_end(sdf, gdf):
    s_beg, s_end = {}, {}
    for sys_id, row in sdf.iterrows():
        s_beg_id = row["sys_beg"]
        s_end_id = row["sys_end"]
        bs, be = gdf.loc[s_beg_id, "gene_beg"], gdf.loc[s_beg_id, "gene_end"]
        es, ee = gdf.loc[s_end_id, "gene_beg"], gdf.loc[s_end_id, "gene_end"]
        assert bs < be, f"Error: {sys_id=} {s_beg_id=} {bs=} {be=}"
        assert es < ee, f"Error: {sys_id=} {s_end_id=} {es=} {ee=}"
        S, E = -1, -1

        if s_beg_id == s_end_id:
            S, E = bs, be
        else:
            if es > be:
                S, E = bs, ee
            elif ee < be:
                print("inverse gene order")
                S, E = es, be
            elif (bs < es) and (be < ee):
                print(
                    f"possibly overlapping genes: {sys_id=} {bs=}, {be=}, {es=}, {ee=}"
                )
                S, E = bs, ee

        assert (
            np.abs(E - S) < 1e5
        ), f"Error: system too long (wrap-around?) | {sys_id=} {S=} {E=}"
        assert S != -1, f"Error unassigned: {sys_id=} {S=} {E=}"
        assert E != -1, f"Error unassigned: {sys_id=} {S=} {E=}"
        # if sys_id in s_beg:
        #     raise ValueError(f"Error: {sys_id=} already in {s_beg=}")
        # if sys_id in s_end:
        #     raise ValueError(f"Error: {sys_id=} already in {s_end=}")

        # iso = row["iso"]
        # if iso == "NZ_CP076693.1":
        #     print(f"{sys_id=} {iso=} {S=} {E=}")
        #     print(f"{s_beg_id=} {s_end_id=}")
        #     print(f"{bs=} {be=} {es=} {ee=}")

        s_beg[sys_id] = S
        s_end[sys_id] = E
    return s_beg, s_end


if __name__ == "__main__":

    args = parse_args()

    gdf = pd.read_csv(args.input_genes, sep="\t")
    gdf.set_index("hit_id", inplace=True)

    sdf = pd.read_csv(args.input_systems, sep="\t", index_col=0)
    sdf["iso"] = sdf["sys_beg"].str.split("_").apply(lambda x: "_".join(x[:-1]))

    s_beg, s_end = system_beg_end(sdf, gdf)
    sdf["sys_beg_bp"] = pd.Series(s_beg)
    sdf["sys_end_bp"] = pd.Series(s_end)

    cols = {"sys_id": "id", "iso": "iso", "sys_beg_bp": "beg", "sys_end_bp": "end"}
    df = sdf.reset_index()[list(cols.keys())].rename(columns=cols)
    df["type"] = "defense_system"
    df.to_csv(args.output_df, index=False)
