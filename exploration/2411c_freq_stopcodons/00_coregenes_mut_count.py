# %%
import json
import pandas as pd
import numpy as np
from Bio import AlignIO, Seq
import gzip
import pathlib
import utils as ut

res_fld = pathlib.Path("res/f00")
res_fld.mkdir(exist_ok=True)


# load gene clusters info
gene_cluster_json = "../../data/panX/data/ST131_ABC/vis/geneCluster.json"
gdf = ut.load_genecluster_json(gene_cluster_json)

# %%
N = 222
core_mask = gdf["count"] == 222
core_mask &= gdf["dupli"] == "no"
gdf["core"] = core_mask
gdf[core_mask]

# 2653 core gene clusters

# %%
# parse alignments


mut_count = {}
mut_possible = {}
aln_stats = {}
for idx, row in gdf[core_mask].iterrows():
    gene_id = idx
    msa_idx = row["msa"]
    aln_file = f"../../data/panX/data/ST131_ABC/vis/geneCluster/{msa_idx}_na_aln.fa.gz"
    with gzip.open(aln_file, "rt") as f:
        aln = AlignIO.read(f, "fasta")
    res = ut.process_alignment(aln)
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
