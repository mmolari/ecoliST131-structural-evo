# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pypangraph as pp
import utils as ut
from collections import defaultdict


len_thr = 500
dset = "ST131_full"
edge_file = f"../../results/{dset}/backbone_joints/asm20-100-5/core_edges.csv"
pangraph_file = f"../../results/{dset}/pangraph/asm20-100-5-polished.json"

# %%
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()
paths = ut.pangraph_to_path_dict(pan)
edges_df = pd.read_csv(edge_file)

# %%


def is_core(node_id):
    return (bdf.loc[node_id, "len"] >= len_thr) and bdf.loc[node_id, "core"]


# core-junctions dataframe
jdf = []
for iso, path in paths.items():
    junctions = ut.path_junction_split(path, is_core)
    for J in junctions:
        edge = J.flanking_edge()
        L = sum(bdf.loc[node.id, "len"] for node in J.center.nodes)
        jdf.append({"iso": iso, "edge": edge.to_unique_str_id(), "len": L})
jdf = pd.DataFrame(jdf)
jdf = jdf.pivot_table(index="iso", columns="edge", values="len")
jdf
# %%

# verify compatibility with edge file
old_edge_df = pd.read_csv(edge_file)
old_core_edges = old_edge_df[old_edge_df["count"] == len(paths)]["edge"].values
old_core_edges = set(
    [ut.Edge.from_str_id(e).to_unique_str_id() for e in old_core_edges]
)
new_core_edges = set(jdf.columns[jdf.notna().sum() == len(paths)])
assert old_core_edges == new_core_edges

# %%

# select core-edges
core_mask = jdf.notna().sum() == len(paths)
core_edges = jdf.columns[core_mask]
jdf_core = jdf[core_edges]

# select accessory-edges
acc_edges = jdf.columns[~core_mask]
jdf_acc = jdf[acc_edges]


# relative genome size
core_genome = jdf_core.sum().sum() / len(paths)
acc_genome = jdf_acc.sum().sum() / len(paths)
print(f"accessory genome in backbone: {core_genome/1e6:.2f} Mbp")
print(f"accessory genome out of backbone: {acc_genome/1e6:.2f} Mbp")
print(f"n. backbone edges {len(core_edges)}")
print(f"n. off-backbone edges {len(acc_edges)}")

mean_nonempty_len = jdf_core[jdf_core > 0].mean()
bins = np.logspace(1, 6, 100)
h = plt.hist(
    mean_nonempty_len,
    log=True,
    bins=bins,
)
plt.xscale("log")
plt.show()
# %%
mask = mean_nonempty_len < 830
mask &= mean_nonempty_len > 739
mean_nonempty_len[mask].value_counts()
# %%
mask = mean_nonempty_len == 766
suspicious = mean_nonempty_len[mask].index.to_numpy()
# %%

sus = []
for iso, path in paths.items():
    junctions = ut.path_junction_split(path, is_core)
    for J in junctions:
        edge = J.flanking_edge()
        if not edge.to_unique_str_id() in suspicious:
            continue
        path = J.center
        if edge.left.id > edge.right.id:
            path = path.invert()
        L_nodes = [bdf.loc[node.id, "len"] for node in path.nodes]
        L = sum(L_nodes)
        if L != 766:
            continue
        N_nodes = len(path.nodes)
        sus.append(
            {
                "iso": iso,
                "edge": edge.to_unique_str_id(),
                "len": L,
                "nodes": N_nodes,
                "L_nodes": L_nodes,
                "path": path,
            }
        )
sus = pd.DataFrame(sus)
sus = sus.sort_values(["edge", "iso"]).reset_index(drop=True)
sus.to_csv("data/suspicious.csv")
sus
# %%
# it's in many genomes
sus["iso"].value_counts()
# %%

edge = "BDLHMQPIBW_f__NKOVRJRZTL_f"
iso_in = "NZ_CP124372.1"
iso_out = "NZ_CP124339.1"

edge = ut.Edge.from_str_id(edge)
bl = pan.blocks[edge.left.id]
br = pan.blocks[edge.right.id]
pos_l_in = [
    (occ[0], occ[2], *pos) for occ, pos in bl.alignment.pos.items() if occ[0] == iso_in
][0]
pos_l_out = [
    (occ[0], occ[2], *pos) for occ, pos in bl.alignment.pos.items() if occ[0] == iso_out
][0]
pos_r_in = [
    (occ[0], occ[2], *pos) for occ, pos in br.alignment.pos.items() if occ[0] == iso_in
][0]
pos_r_out = [
    (occ[0], occ[2], *pos) for occ, pos in br.alignment.pos.items() if occ[0] == iso_out
][0]
print(pos_l_in)
print(pos_r_in)
print(pos_l_out)
print(pos_r_out)

# extract sequence from iso_in
from Bio import SeqIO

with open(f"../../data/fa/{iso_in}.fa") as f:
    genome_in = SeqIO.read(f, "fasta")
if pos_l_in[3] < pos_r_in[2]:
    seq_sus = genome_in.seq[pos_l_in[3] : pos_r_in[2]]
else:
    seq_sus = genome_in.seq[pos_r_in[3] : pos_l_in[2]]
print(seq_sus)

# %%

import subprocess


def disentangle_pos(pos_l, pos_r):
    sdf = defaultdict(dict)
    for occ, pos in pos_l.items():
        sdf[occ[0]]["lb"] = pos[0]
        sdf[occ[0]]["le"] = pos[1]
        sdf[occ[0]]["ls"] = occ[2]
    for occ, pos in pos_r.items():
        sdf[occ[0]]["rb"] = pos[0]
        sdf[occ[0]]["re"] = pos[1]
        sdf[occ[0]]["rs"] = occ[2]
    return pd.DataFrame(sdf).T


def load_genome(iso):
    with open(f"../../data/fa/{iso}.fa") as f:
        genome = SeqIO.read(f, "fasta")
    return genome.seq


def extract_seq(seq, delta, pos):
    avg_l = (pos["lb"] + pos["le"]) / 2
    avg_r = (pos["rb"] + pos["re"]) / 2
    assert np.abs(avg_l - avg_r) < 1e6, f"weird positions: {pos}"
    fwd = avg_l < avg_r
    beg = pos["le"] if fwd else pos["re"]
    end = pos["rb"] if fwd else pos["lb"]
    seq = seq[beg - delta : end + delta]
    if not fwd:
        seq = seq.reverse_complement()
    return seq


delta = 500
for edge, sdf in sus.groupby("edge"):
    print(sdf)

    edge = ut.Edge.from_str_id(edge)
    pos_l = pan.blocks[edge.left.id].alignment.pos
    pos_r = pan.blocks[edge.right.id].alignment.pos
    pos = disentangle_pos(pos_l, pos_r)
    records = []
    for iso, row in pos.iterrows():
        genome = load_genome(iso)
        seq = extract_seq(genome, delta, row)
        records.append(SeqIO.SeqRecord(seq, id=iso))
    fname = f"data/IS/{edge.to_unique_str_id()}.fa"
    SeqIO.write(records, fname, "fasta")

    # align usign mafft
    fname_out = f"data/IS/{edge.to_unique_str_id()}.aln.fa"
    subprocess.run(
        f"mafft --quiet --auto {fname} > {fname_out}",
        shell=True,
    )
    # remove unaligned file
    subprocess.run(f"rm {fname}", shell=True)
    break
# %%
sus_isos = sus["iso"].unique()
tree_file = f"../../results/{dset}/pangraph/asm20-100-5-filtered-coretree.nwk"
from Bio import Phylo

tree = Phylo.read(tree_file, "newick")
tree.root_at_midpoint()
tree.ladderize()

# %%
fig, ax = plt.subplots(figsize=(10, 20))
Phylo.draw(
    tree,
    label_func=lambda x: x.name if x.name in sus_isos else "",
    do_show=False,
    axes=ax,
)
plt.tight_layout()
# plt.savefig("figs/missed_accessory/tree.png")
plt.show()

# %%
dset = "ST131_full"
tree_file = f"../../results/{dset}/pangraph/asm20-100-5-filtered-coretree.nwk"
from Bio import Phylo

tree = Phylo.read(tree_file, "newick")
tree.root_at_midpoint()
tree.ladderize()
# %%
import matplotlib.pyplot as plt

names = [l.name for l in tree.get_terminals()]
fig, ax = plt.subplots(figsize=(10, 30))
Phylo.draw(
    tree,
    label_func=lambda x: x.name if x.name in names else "",
    do_show=False,
    axes=ax,
)
plt.tight_layout()
# plt.savefig("figs/missed_accessory/tree.png")
plt.show()
# %%
