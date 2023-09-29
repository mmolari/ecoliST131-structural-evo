import pandas as pd
import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--blast_tsv",
        type=str,
        help="input blast file",
    )
    parser.add_argument(
        "--tsv_out",
        type=str,
        help="output file",
    )
    parser.add_argument(
        "--min_cov",
        type=float,
        help="minimal coverage to assign an allele",
    )
    parser.add_argument(
        "--min_id",
        type=float,
        help="minimal identity to assign an allele",
    )
    return parser.parse_args()


class match:
    def __init__(self) -> None:
        self.locus = ""
        self.allele = ""
        self.match_type = ""
        self.cov = ""
        self.sim = ""
        self.matches = ""
        self.aln_L = ""
        self.iso = ""
        self.allele_L = ""

    def assign(self, row) -> None:
        self.allele = row["ref_id"].split("_")[1]
        self.cov = row["cov"]
        self.sim = row["sim"]
        self.matches = row["matches"]
        self.aln_L = row["aln_L"]
        self.allele_L = row["rL"]

    def to_tsv_string(self) -> str:
        return "\t".join(
            [
                self.iso,
                self.locus,
                self.match_type,
                self.allele,
                str(self.cov),
                str(self.sim),
                str(self.matches),
                str(self.aln_L),
                str(self.allele_L),
            ]
        )


def import_df(fname):
    df = pd.read_csv(fname, sep="\t", header=None)
    if len(df) == 0:
        return None
    df.columns = [
        "ref_id",
        "rL",
        "aln_L",
        "matches",
        "qry_id",
        "qs",
        "qe",
        "qstrand",
    ]
    df.sort_values("matches", ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


def find_match(df, m, min_cov, min_id):
    df["cov"] = df["aln_L"] / df["rL"]
    df["sim"] = df["matches"] / df["aln_L"]

    # above minimal coverage
    mask = df["cov"] >= min_cov
    sdf = df[mask]

    # above minimal identity
    mask = sdf["sim"] >= min_id
    sdf = sdf[mask]

    # find best match
    sdf.sort_values("matches", ascending=False, inplace=True)
    sdf.reset_index(drop=True, inplace=True)

    if len(sdf) == 0:
        m.match_type = "-"
        return m
    elif len(sdf) == 1:
        m.match_type = "1"
        m.assign(sdf.iloc[0])
        return m

    # find two best matches
    m1 = sdf.iloc[0]
    m2 = sdf.iloc[1]
    m.assign(m1)

    identical = m1["matches"] == m1["rL"]
    best = m1["matches"] > m2["matches"]
    equal = m1["matches"] == m2["matches"]
    if identical:
        m.match_type = "+"
    elif best:
        m.match_type = "~"
    elif equal:
        m.match_type = "?"
    else:
        raise ValueError("wrong match order")

    return m


if __name__ == "__main__":
    args = parse_args()
    iso = args.blast_tsv.split("/")[-1].split(".")[0]
    locus = args.blast_tsv.split("/")[-2]
    m = match()
    m.iso = iso
    m.locus = locus
    m.match_type = "x"

    fsize = os.path.getsize(args.blast_tsv)
    if fsize > 0:
        df = import_df(args.blast_tsv)
        if df is not None:
            m = find_match(df, m, args.min_cov, args.min_id)

    with open(args.tsv_out, "w") as f:
        f.write(m.to_tsv_string())
        f.write("\n")
