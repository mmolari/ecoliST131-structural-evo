# %%
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hochhauser", type=str)
    parser.add_argument("--gbk", type=str, nargs="+")
    parser.add_argument("--out", type=str)
    return parser.parse_args()


def load_hochhauser_genes(fname):
    hdf = pd.read_csv(fname, index_col=0)
    cols = [
        "Gene symbol of upstream gene",
        "Gene symbol of downstream gene",
    ]
    genes = set(hdf[cols].values.flatten())
    return genes


def find_genes(genes, fname):
    gene_ctr = defaultdict(int)
    res = []
    with open(fname, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            iso = record.id
            assert fname.endswith(iso + ".gbk")
            for feature in record.features:
                if feature.type == "CDS":
                    if "gene" in feature.qualifiers:
                        gene = feature.qualifiers["gene"][0]
                        if gene in genes:
                            gene_ctr[gene] += 1
                            res.append(
                                {
                                    "id": f"{gene}|{iso}|{gene_ctr[gene]}",
                                    "iso": iso,
                                    "beg": feature.location.start,
                                    "end": feature.location.end,
                                    "strand": feature.location.strand,
                                    "type": gene,
                                }
                            )
    return res


if __name__ == "__main__":

    args = parse_args()
    genes = load_hochhauser_genes(args.hochhauser)
    gdf = []
    for fname in args.gbk:
        gdf += find_genes(genes, fname)
    gdf = pd.DataFrame(gdf)
    gdf.to_csv(args.out, index=False)

# %%
