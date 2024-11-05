# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
from Bio import Phylo
from utils import mugration_inference, name_tree_nodes
import json
import subprocess
from multiprocessing import Pool


res_fld = pathlib.Path("res")

# load gene counts
GC = pd.read_csv(res_fld / "gene_counts.csv", index_col=0)


# load tree:
pangraph_tree_file = (
    "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
)
tree = Phylo.read(pangraph_tree_file, "newick")
name_tree_nodes(tree)
strains = [l.name for l in tree.get_terminals()]

# load gene clusters info
gene_cluster_json = "../../data/panX/data/ST131_ABC/vis/geneCluster.json"

with open(gene_cluster_json, "r") as f:
    gdf = pd.DataFrame(json.load(f))
gdf["divers"] = gdf["divers"].astype(float)
gdf["count"] = gdf["count"].astype(int)
gdf["geneLen"] = gdf["geneLen"].astype(int)
gdf["event"] = gdf["event"].astype(int)

# %%
mask_core = (GC == GC.iloc[0]).all(axis=0)
GC_acc = GC.loc[:, ~mask_core]


chunk_size = 10
chunk_idxs = list(range(0, GC_acc.shape[1], chunk_size)) + [GC_acc.shape[1]]
chunk_idxs = np.unique(chunk_idxs)  # remove duplicates


def process_chunk_i(i):
    beg, end = chunk_idxs[i], chunk_idxs[i + 1]
    chunk = GC_acc.iloc[:, beg:end].copy()
    # convert numbers to letters
    chunk = chunk.map(lambda x: chr(x + 65))
    mugr_res = mugration_inference(tree, chunk)
    return mugr_res


with Pool(14) as p:
    res = p.map(process_chunk_i, range(len(chunk_idxs) - 1))

# join dictionaries
mugr_res = {}
for r in res:
    mugr_res.update(r)

# save gzipped json file
with open(res_fld / "mugr_res.json", "w") as f:
    json.dump(mugr_res, f)
# gzip
subprocess.run(["gzip", res_fld / "mugr_res.json"])

# %%


def parse_events(mugr_res):
    E = []
    for k, v in mugr_res.items():
        ev = v["pa_pattern"]["events"]
        for branch, event in ev:
            after, before = event.split("|")
            delta = ord(after) - ord(before)
            E.append((k, branch, delta))

    E = pd.DataFrame(E, columns=["gcl", "branch", "delta"])
    return E


# extract events
E = parse_events(mugr_res)
E["type"] = E["delta"].map(lambda x: "gain" if x > 0 else "loss")
E.sort_values("gcl", inplace=True)

# save events
E.to_csv(res_fld / "mugr_events.csv", index=False)

# %%
