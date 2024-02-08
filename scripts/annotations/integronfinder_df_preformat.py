import pandas as pd
import numpy as np
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_df", required=True)
    parser.add_argument("--output_df", required=True)
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()
    df = pd.read_csv(args.input_df, sep="\t")

    df["integron_n"] = df["ID_integron"].str.replace("integron_", "").astype(int)

    beg = df.groupby(["ID_replicon", "integron_n"])["pos_beg"].min()
    end = df.groupby(["ID_replicon", "integron_n"])["pos_end"].max()
    tp = df.groupby(["ID_replicon", "integron_n"])["type"].first()

    idf = pd.DataFrame({"beg": beg, "end": end, "type": tp})

    for idx, row in idf.iterrows():
        if np.abs(row["beg"] - row["end"]) > 1e6:
            # possible error with beg/end definition.
            # if wrapping around the genome then it must be defined in a
            # different way
            raise ValueError(f"Integron {idx} is too big")
    idf["len"] = idf["end"] - idf["beg"] + 1

    idf = idf.reset_index()
    idf.rename(columns={"ID_replicon": "iso"}, inplace=True)
    idf["id"] = idf["iso"] + "|" + idf["integron_n"].astype(str) + "|" + idf["type"]
    cols = ["id", "iso", "beg", "end", "type"]
    idf = idf[cols]
    idf.set_index("id", inplace=True, verify_integrity=True)
    idf.to_csv(args.output_df)
