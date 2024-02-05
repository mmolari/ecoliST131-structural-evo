import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_df", required=True)
    parser.add_argument("--output_df", required=True)
    return parser.parse_args()


# cols = [
#     "seqID",
#     "family",
#     "cluster",
#     "isBegin",
#     "isEnd",
#     "isLen",
#     "ncopy4is",
#     "start1",
#     "end1",
#     "start2",
#     "end2",
#     "score",
#     "irId",
#     "irLen",
#     "nGaps",
#     "orfBegin",
#     "orfEnd",
#     "strand",
#     "orfLen",
#     "E-value",
#     "E-value4copy",
#     "type",
#     "ov",
#     "tir",
# ]
# columns_final = ["id", "type", "iso", "beg", "end", "strand"]


if __name__ == "__main__":

    args = parse_args()
    df = pd.read_csv(args.input_df, sep="\t")

    cols = {
        "seqID": "iso",
        "family": "id",
        "isBegin": "beg",
        "isEnd": "end",
        "strand": "strand",
    }

    sdf = df[list(cols.keys())].rename(columns=cols)
    sdf["type"] = "IS"
    sdf["id"] = sdf["iso"] + "|" + sdf["id"] + "|" + sdf.index.astype(str)
    sdf.set_index("id", inplace=True, verify_integrity=True)

    sdf.to_csv(args.output_df)
