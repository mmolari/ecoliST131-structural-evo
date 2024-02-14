import pandas as pd
from Bio import SeqIO
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_gene_df", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--output_gene_df", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    genes_df_file = args.input_gene_df
    prot_file = args.proteins
    out_file = args.output_gene_df

    df = pd.read_csv(genes_df_file, sep="\t")
    with open(prot_file, "r") as f:
        prots = list(SeqIO.parse(f, "fasta"))
    prots = {p.id: p.description for p in prots}

    df["gene_beg"] = -1
    df["gene_end"] = -1
    df["gene_strand"] = ""
    for i, row in df.iterrows():
        idx = row["hit_id"]
        descr = prots[idx]
        beg, end, strand = [int(descr.split("#")[x].strip()) for x in [1, 2, 3]]
        df.at[i, "gene_beg"] = beg
        df.at[i, "gene_end"] = end
        df.at[i, "gene_strand"] = strand
    df.to_csv(out_file, sep="\t", index=False)
