# %%
import json
import pandas as pd
import numpy as np
from Bio import AlignIO, Seq
import pathlib
import utils as ut

res_fld = pathlib.Path("res/f02")
res_fld.mkdir(exist_ok=True)

load_fld = pathlib.Path("res/f01")

gene_cluster_json = "../../data/panX/data/ST131_ABC/vis/geneCluster.json"
gc = ut.load_genecluster_json(gene_cluster_json)
# %%
gdf = pd.read_csv(load_fld / "stats.csv", index_col=0)
gdf.set_index("gid", inplace=True)
gdf
# %%


mut_count = {}
mut_possible = {}
aln_stats = {}
for idx, row in gdf.iterrows():
    print(idx)
    gene_id = idx
    aln_file = load_fld / f"{gene_id}_aln.fa"
    with open(aln_file, "r") as f:
        aln = AlignIO.read(f, "fasta")
    A = np.array(aln)

    # remove columns with majority gaps
    mask = (A == "-").sum(axis=0) < 0.5 * len(A)
    A = A[:, mask]

    # gene length
    gene_L = gc.loc[gene_id, "geneLen"]
    if gene_L % 3 == 0:
        A = A[:, : gene_L - 3]
    else:
        A = A[:, :-3]

    res = ut.process_alignment(A)
    mut_count[gene_id] = res["mut_count"]
    mut_possible[gene_id] = res["mut_possible"]
    aln_stats[gene_id] = res["aln_stats"]

mut_count = pd.DataFrame(mut_count).fillna(0).astype(int)
mut_possible = pd.DataFrame(mut_possible).fillna(0).astype(int)
aln_stats = pd.DataFrame(aln_stats).fillna(0).astype(int)

# %%

mut_count.to_csv(res_fld / "mut_count.csv")
mut_possible.to_csv(res_fld / "mut_possible.csv")
aln_stats.to_csv(res_fld / "aln_stats.csv")

# %%
