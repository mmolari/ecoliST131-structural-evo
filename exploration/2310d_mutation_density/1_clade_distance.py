# %%
from Bio import AlignIO, SeqIO, Phylo
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

dset = "ST131_full"
subfld = f"../../results/{dset}/pangraph"
tree_file = f"{subfld}/asm20-100-5-filtered-coretree.nwk"
alignment_file = f"{subfld}/asm20-100-5-alignment/corealignment.fa"
filtered_alignment_file = f"{subfld}/asm20-100-5-alignment/filtered_corealignment.fa"
# aln_info = f"{subfld}/asm20-100-5-alignment/corealignment_info.json"
# filt_aln_info = f"{subfld}/asm20-100-5-alignment/filtered_corealignment_info_size.json"


def consensus(A):
    """given an alignment matrix returns the consensus"""

    def site_consensus(l):
        lett, ct = np.unique(l, return_counts=True, axis=0)
        return lett[np.argmax(ct)]

    return np.apply_along_axis(site_consensus, 0, A)


filt = False
filt_label = "filtered" if filt else "unfiltered"
if filt:
    aln = AlignIO.read(filtered_alignment_file, "fasta")
else:
    aln = AlignIO.read(alignment_file, "fasta")
A = np.array(aln)
iso_aln_idx = {iso.id: i for i, iso in enumerate(aln)}

tree = Phylo.read(tree_file, "newick")
tree.root_at_midpoint()
tree.ladderize()

# %%
# assign_clades
clade = {}
for node in tree.root[0].get_terminals():
    clade[node.name] = "A"
for node in tree.root[1].get_terminals():
    clade[node.name] = "B"
for node in tree.root[1][2].get_terminals():
    clade[node.name] = "C"

fig, ax = plt.subplots(1, 1, figsize=(4, 15))
Phylo.draw(
    tree,
    label_func=lambda x: x.name if x.is_terminal() else "",
    do_show=False,
    label_colors=lambda x: {"A": "red", "B": "blue", "C": "green"}[clade[x]]
    if len(x) > 0
    else "black",
    axes=ax,
)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.tight_layout()
plt.savefig("figs/clade_tree.png", facecolor="w")
plt.show()

# %%
W = 50
cons, D = {}, {}
for cl in ["A", "B", "C"]:
    mask = np.array([clade[a.id] == cl for a in aln])
    cons[cl] = consensus(A[mask])
    D[cl] = np.diff((A != cons[cl]).cumsum(axis=1)[:, ::W])
    nW = D[cl].shape[1]

# %%
fig, axs = plt.subplots(
    1, 2, figsize=(8, 60), sharey=True, gridspec_kw={"width_ratios": [1, 3]}
)

ax = axs[0]
Phylo.draw(
    tree,
    label_func=lambda x: clade[x.name] if x.is_terminal() else "",
    do_show=False,
    label_colors={"A": "red", "B": "blue", "C": "green", "": "black"},
    axes=ax,
)

ax = axs[1]

for n, iso in enumerate(tree.get_terminals()):
    y0 = (n + 1) * np.ones(nW)

    for cl in ["A", "B", "C"]:
        idx = iso_aln_idx[iso.name]
        y = -D[cl][idx] / W
        x = np.arange(nW) * W + W / 2
        ax.plot(x, y0 + y, color={"A": "red", "B": "blue", "C": "green"}[cl])
ax.set_title(f"avg distance to clade in {W} bp window")
plt.tight_layout()
plt.savefig(f"figs/clade_distance_{filt_label}_W_{W}.png", dpi=300, facecolor="w")
plt.show()

# %%
# assign similarity to iso:
thr = 0.10
D_assign = {}
for iso in tree.get_terminals():
    idx = iso_aln_idx[iso.name]
    res = defaultdict(set)
    for cl in ["A", "B", "C"]:
        for n, v in enumerate(D[cl][idx] < thr):
            if v:
                res[n * W + W / 2] |= set(cl)
    D_assign[iso.name] = res
D_df = pd.DataFrame(D_assign)
D_df = D_df.applymap(lambda x: "|".join(sorted(x)) if not pd.isna(x) else "o")


colors = {
    "A|B|C": "white",
    "A|B": "magenta",
    "A|C": "yellow",
    "B|C": "cyan",
    "A": "red",
    "B": "blue",
    "C": "green",
    "o": "black",
}

fig, axs = plt.subplots(
    1, 2, figsize=(15, 30), sharey=True, gridspec_kw={"width_ratios": [1, 5]}
)

ax = axs[0]
Phylo.draw(
    tree,
    label_func=lambda x: clade[x.name] if x.is_terminal() else "",
    do_show=False,
    label_colors={"A": "red", "B": "blue", "C": "green", "": "black"},
    axes=ax,
)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

ax = axs[1]
for n, iso in enumerate(tree.get_terminals()):
    idx = iso_aln_idx[iso.name]
    y = n + 1
    Xs = np.sort(D_df.index)
    for x in Xs:
        ax.barh(y, W, left=x - W / 2, color=colors[D_df.loc[x, iso.name]])


plt.tight_layout()
plt.savefig(f"figs/clade_tree_color_{filt_label}_thr_{thr}.png", facecolor="w")
plt.show()

# %%
