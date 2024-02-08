import pandas as pd
from Bio import SeqIO
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastas", nargs="+", help="Fasta file(s) to get length from")
    parser.add_argument("--df", "--output", help="Output file name")
    return parser.parse_args()


def get_info(fa):
    with open(fa) as f:
        record = SeqIO.read(f, "fasta")
    return record.id, len(record.seq)


if __name__ == "__main__":
    args = parse_args()

    data = []
    for fa in args.fastas:
        idx, L = get_info(fa)
        data.append([idx, L])
    df = pd.DataFrame(data, columns=["iso", "length"])
    df.to_csv(args.df, index=False)
