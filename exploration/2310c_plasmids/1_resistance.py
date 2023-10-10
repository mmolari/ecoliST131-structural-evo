import argparse
import pathlib
import json
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="assign plasmid resistance summary to corresponding isolates."
    )
    parser.add_argument("--plasmid_json", type=str, help="plasmid resistance file")
    parser.add_argument(
        "--threshold_id",
        type=float,
        help="threshold identity fraction for resistance genes",
    )
    parser.add_argument("--resistance_tsv", type=str, help="resistance tsv file")
    parser.add_argument("--out_csv", type=str, help="output csv file")
    return parser.parse_args()


def parse_element(x, thr):
    """
    Assign `.` = 0
    and `100.00;100;00;100.00` = 3
    """
    if x == ".":
        return 0
    x = str(x)
    els = x.split(";")
    ct = 0
    for el in els:
        if float(el) > thr * 100.0:
            ct += 1
    return ct


def load_res_df(fname, thr):
    df = pd.read_csv(fname, sep="\t", index_col=0)
    # transform filename to accession number
    df.index = df.index.map(lambda x: pathlib.Path(x).stem)
    df.drop(columns="NUM_FOUND", inplace=True)
    # parse elements
    df = df.applymap(lambda x: parse_element(x, thr))
    return df


def assing_to_chromosome(df, plasmids_dict):
    """
    Assign plasmid resistance to chromosome
    """
    rows = []
    for acc, plasmids in plasmids_dict.items():
        row = df.loc[plasmids].sum(axis=0)
        row.name = acc
        rows.append(row)
    df_chr = pd.DataFrame(rows)
    return df_chr


if __name__ == "__main__":
    args = parse_args()
    # load resistance df
    df = load_res_df(args.resistance_tsv, args.threshold_id)
    # load plasmid dictionary
    plasmids_dict = json.load(open(args.plasmid_json))
    # assign plasmid resistance to isolates
    df_chr = assing_to_chromosome(df, plasmids_dict)
    # save
    df_chr.to_csv(args.out_csv)
