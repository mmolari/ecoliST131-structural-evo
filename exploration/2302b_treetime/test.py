# %%

import json
import treetime

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sp
import seaborn as sns


from Bio import Phylo
from collections import Counter


tree_file = "../../results/ST131/pangraph/asm20-100-5-coretree.nwk"
aln_file = "../../results/ST131/pangraph/asm20-100-5-alignment/corealignment.fa"
aln_info_file = (
    "../../results/ST131/pangraph/asm20-100-5-alignment/corealignment_info.json"
)

with open(aln_info_file, "r") as f:
    info = json.load(f)
N_snps = info["n. snps"]
aln_L = info["n. consensus"] + info["n. snps"]
factor = info["n. snps"] / aln_L

# %%

tree = Phylo.read(tree_file, format="newick")
Phylo.draw(tree.root, label_func=lambda x: "", show_confidence=False)


# instantiate treetime
myTree = treetime.TreeAnc(
    gtr="Jukes-Cantor", tree=tree, aln=aln_file, verbose=0, seq_len=aln_L
)
myTree.tree.root.branch_length = 0.0


# %%
myTree.infer_ancestral_sequences(prune_short=True)


# %%
mut_count = Counter()
bl = []
nm = []
for node in myTree.tree.find_clades():
    muts = node.mutations
    mut_count.update([m[1] for m in muts])
    # print(node, node.branch_length, len(muts))
    bl.append(node.branch_length)
    nm.append(len(muts))
    # print(muts)

slope = sp.stats.linregress(bl, nm).slope

plt.figure(figsize=(3, 3))
plt.scatter(bl, nm, alpha=0.1)
x = np.unique(np.sort(bl))
plt.plot(x, x * slope, "k--")
plt.xscale("symlog", linthresh=1e-7)
plt.yscale("symlog", linthresh=1)
plt.xlabel("branch length")
plt.ylabel("n. mutations")
plt.tight_layout()
plt.show()
# %%

plt.figure(figsize=(8, 1.5))
x = np.sort(list(mut_count.keys()))
y = np.array([mut_count[k] for k in x])
plt.scatter(x, y > 1, alpha=0.2, marker="|")
# plt.yscale("symlog")
plt.xlabel("core genome restricted alignment (bp)")
plt.ylabel("site mutations")
plt.yticks([0, 1], ["single", "multiple"])
plt.tight_layout()
plt.show()

# %%


def partition_entropy(ps, L):
    x = np.sort(ps + [0, L])
    deltas = np.diff(x) / L
    return sp.stats.entropy(deltas, base=np.exp(1))


df_supporting = []
nodes = list(myTree.tree.find_clades())
N_nodes = len(nodes)

bins = np.arange(0, N_snps + 50, 25)

fig, axs = plt.subplots(N_nodes, 1, sharex=True, figsize=(8, N_nodes * 1.0))
# fig, axs = plt.subplots(1, 1, sharex=True, figsize=(8, 2.0))
# axs = [axs]

for idx, node in enumerate(sorted(nodes, key=lambda x: x.branch_length, reverse=True)):
    muts = node.mutations
    ps = [m[1] for m in muts]
    pinc = [p for p in ps if mut_count[p] > 1]
    psupp = [p for p in ps if mut_count[p] == 1]
    res = {
        "node": node.name,
        "len": node.branch_length,
        "n. supporting": len(psupp),
        "n. double": len(pinc),
        "part. entropy": partition_entropy(psupp, N_snps),
    }
    print(res)
    df_supporting.append(res)
    # if len(ps) == 0:
    #     continue
    axs[idx].hist(psupp, bins=bins, cumulative=False)
    axs[idx].hist(pinc, bins=bins, cumulative=False, weights=[-1] * len(pinc))
    axs[idx].set_title(
        f"node: {node.name} -- len: {node.branch_length} -- n. muts {len(psupp)}"
    )
    axs[idx].grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("figs/branch_support_vs_aln.png")
plt.show()
df_supporting = pd.DataFrame(df_supporting)

# %%


plt.figure(figsize=(3, 3))
plt.scatter(df_supporting.len, df_supporting["n. supporting"], alpha=0.2)
x = np.unique(np.sort(bl))
plt.plot(x, x * slope, "k--", alpha=0.3)
plt.xscale("symlog", linthresh=1e-7)
plt.yscale("symlog", linthresh=1)
plt.xlabel("branch length")
plt.ylabel("n. single mutations")
plt.tight_layout()
plt.show()

# %%

df = []
for n in [1, 2, 3, 4, 5, 10, 15, 20, 30, 40, 50, 100, 500, 1000, 5000, 10000]:
    for _ in range(100):
        x = np.random.rand(n)
        x = np.sort(list(x) + [0, 1])
        x = np.diff(x)
        df.append({"n": n, "H": sp.stats.entropy(x, base=np.exp(1))})
df = pd.DataFrame(df)
df["dH"] = np.log(df["n"] + 1) - df["H"]

gamma = 0.577215
sns.catplot(df, x="n", y="dH", alpha=0.2)
plt.axhline(1 - gamma, color="k")
plt.show()
# %%

sdf = df_supporting[df_supporting["n. supporting"] > 0].copy()
sdf["n"] = sdf["n. supporting"]
sdf = sdf.rename(columns={"part. entropy": "H"})
sdf["dH"] = np.log(sdf["n"] + 1) - sdf["H"]

ns = np.sort(sdf["n"].unique())
df = []
for n in ns:
    for _ in range(100):
        x = np.random.rand(n)
        x = np.sort(list(x) + [0, 1])
        x = np.diff(x)
        df.append({"n": n, "H": sp.stats.entropy(x, base=np.exp(1))})
df = pd.DataFrame(df)
df["dH"] = np.log(df["n"] + 1) - df["H"]

# %%

df["type"] = "rand"
sdf["type"] = "data"
fdf = pd.concat([df, sdf], axis=0)

import seaborn as sns


fig, ax = plt.subplots(1, 1, figsize=(10, 5))
plt.axhline(1 - gamma, color="k")
sns.stripplot(fdf, x="n", y="dH", hue="type", ax=ax)
plt.xticks(rotation=90)
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(10, 5))
# plt.axhline(1 - gamma, color="k")
sns.stripplot(fdf, x="n", y="H", hue="type", ax=ax)
plt.xticks(rotation=90)
plt.show()

# %%
fdf["n_inf"] = np.exp(fdf["H"] + 1 - gamma) - 1
# %%
sns.scatterplot(fdf, x="n", y="n_inf", hue="type")
plt.xscale("log")
plt.yscale("log")
plt.grid(which="major", alpha=0.3)
plt.grid(which="minor", alpha=0.1)
plt.show()
# %%
mask = fdf["dH"] > 1.5
nodes_df = fdf[mask]

# %%
fig, ax = plt.subplots(1, 1, figsize=(5, 10))
Phylo.draw(
    myTree.tree,
    label_func=lambda x: x.name.removeprefix("NODE_00000").removesuffix(".00") if x.name in nodes_df["node"].to_list() else "",
    label_colors=lambda x: "r",
    axes=ax
)
# %%
