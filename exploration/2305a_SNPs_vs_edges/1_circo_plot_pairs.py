# %%

import utils as ut

import numpy as np
from Bio import SeqIO, AlignIO, Phylo
from pycirclize import Circos

import pathlib

import matplotlib as mpl
import matplotlib.pyplot as plt

fig_p = pathlib.Path("fig")
fig_p.mkdir(exist_ok=True)


def svfig(svname):
    for k in ["pdf", "png"]:
        plt.savefig(fig_p / f"{svname}.{k}", dpi=300, facecolor="w")


prefix = "../../results/ST131/pangraph"

aln_file = f"{prefix}/asm20-100-5-alignment/filtered_corealignment.fa"
tree_file = f"{prefix}/asm20-100-5-filtered-coretree.nwk"
pangraph_file = f"{prefix}/asm20-100-5-polished.json"

tree = Phylo.read(tree_file, format="newick")
tree.ladderize()
str_order = [l.name for l in tree.get_terminals()]

aln = AlignIO.read(aln_file, "fasta")
A = np.array(aln)
aln_strains = [s.name for s in aln]
# %%

# extract SNPs pairs
snps_pairs = ut.aln_to_branchings(aln_file)
snps_pairs = [b for b in snps_pairs if len(b) == 2]
# %%
block_pa = ut.extract_block_PA_df(pangraph_file, len_thr=100, remove_dupl=False)

# %%

bl_pairs = []
for bl, pa in block_pa.iteritems():
    vc = pa.value_counts()
    if len(vc) != 2:
        continue
    if vc.min() != 2:
        continue
    idxs = pa.index[pa == vc.index[1]]
    bl_pairs.append(list(idxs))
# %%

# extract edge dictionary
E_dict = ut.extract_edge_dictionary(pangraph_file, len_thr=100, remove_dupl=False)
strains_set = set(aln_strains)
N = len(strains_set)
# %%
edge_pairs = []
for k, s in E_dict.items():
    if len(s) == 2:
        # check that the block is present in both
        # b1, b2 = k.b[0].id, k.b[1].id
        # if block_pa[[b1, b2]].loc[s].all().any():
        edge_pairs.append(list(s))
    # if len(s) == N - 2:
    # edge_pairs.append(list(strains_set - s))


# %%

# Initialize circos sector with tree size
circos = Circos(sectors={"Tree": tree.count_terminals()})
sector = circos.sectors[0]

# Plot tree with inner style
track = sector.add_track((75, 100))
track.axis(ec="lightgrey")
track.tree(
    tree,
    outer=False,
    use_branch_length=True,
    leaf_label_size=0,
    line_kws=dict(color="black", lw=1),
    text_kws=dict(color="gray"),
)

circos.text(
    "\n".join(
        [
            f"n. pairwise",
            f"snps={len(snps_pairs)}",
            f"blocks={len(bl_pairs)}",
            f"edges={len(edge_pairs)}",
        ]
    ),
    size=10,
    deg=315,
    r=120,
)

for p1, p2 in snps_pairs:
    i1 = str_order.index(p1) + 0.5
    i2 = str_order.index(p2) + 0.5
    reg1 = ("Tree", i1, i1)
    reg2 = ("Tree", i2, i2)
    circos.link(reg1, reg2, color="red", lw=0.5, alpha=0.2, height_ratio=0.1)

for p1, p2 in bl_pairs:
    i1 = str_order.index(p1) + 0.5
    i2 = str_order.index(p2) + 0.5
    reg1 = ("Tree", i1, i1)
    reg2 = ("Tree", i2, i2)
    circos.link(reg1, reg2, color="blue", lw=0.5, alpha=0.05, height_ratio=0.2)

for p1, p2 in edge_pairs:
    i1 = str_order.index(p1) + 0.5
    i2 = str_order.index(p2) + 0.5
    reg1 = ("Tree", i1, i1)
    reg2 = ("Tree", i2, i2)
    circos.link(reg1, reg2, color="green", lw=0.5, alpha=0.05, height_ratio=0.23)


fig = circos.plotfig()
plt.tight_layout()
svfig("pairs_on_tree")
plt.show()
# %%
N = len(strains_set)
Ms = np.zeros((N, N))
Mb = np.zeros((N, N))
Me = np.zeros((N, N))
for p1, p2 in snps_pairs:
    i1 = str_order.index(p1)
    i2 = str_order.index(p2)
    Ms[i1, i2] += 1
    Ms[i2, i1] += 1

for p1, p2 in bl_pairs:
    i1 = str_order.index(p1)
    i2 = str_order.index(p2)
    Mb[i1, i2] += 1
    Mb[i2, i1] += 1

for p1, p2 in edge_pairs:
    i1 = str_order.index(p1)
    i2 = str_order.index(p2)
    Me[i1, i2] += 1
    Me[i2, i1] += 1

# for M in [Ms, Mb, Me]:
#     M[:, :] = np.log10(M)

# %%

fig, axs = plt.subplots(1, 3, figsize=(12, 3.6), sharex=True, sharey=True)

for ax, M, lab in zip(axs, [Ms, Mb, Me], ["SNPs", "Blocks", "Edges"]):
    norm = mpl.colors.LogNorm(vmin=1, vmax=M.max())

    g = ax.matshow(M, cmap="cool", norm=norm)
    plt.colorbar(g, ax=ax, shrink=0.8)
    ax.set_title(lab)

    ax.set_xlim(-1, N)
    ax.set_ylim(N, -1)

    # set major ticks
    ax.set_xticks(np.arange(0, N, 10))
    ax.set_yticks(np.arange(0, N, 10))

    # set minor ticks
    ax.set_xticks(np.arange(0, N, 5), minor=True)
    ax.set_yticks(np.arange(0, N, 5), minor=True)
    ax.grid(which="minor", alpha=0.15)
    ax.grid(alpha=0.4)


plt.tight_layout()
svfig("pairs_on_matrices")
plt.show()


# %%

Msb = np.copy(Ms)
Mse = np.copy(Ms)
Mbe = np.copy(Mb)
# make the upper triangular part equal to Mb
Msb[np.triu_indices(N, 1)] = Mb[np.triu_indices(N, 1)]
Mse[np.triu_indices(N, 1)] = Me[np.triu_indices(N, 1)]
Mbe[np.triu_indices(N, 1)] = Me[np.triu_indices(N, 1)]

# display the

# %%
fig, axs = plt.subplots(1, 3, figsize=(12, 3.6), sharex=True, sharey=True)

for ax, M, lab in zip(
    axs, [Msb, Mse, Mbe], [["SNPs", "blocks"], ["SNPs", "edges"], ["blocks", "edges"]]
):
    norm = mpl.colors.LogNorm(vmin=1, vmax=M.max())
    g = ax.matshow(M, cmap="cool", norm=norm)
    plt.colorbar(g, ax=ax, shrink=0.8)
    ax.set_ylabel(lab[0])
    ax.set_title(lab[1])

    ax.set_xlim(-1, N)
    ax.set_ylim(N, -1)

    # set major ticks
    ax.set_xticks(np.arange(0, N, 10))
    ax.set_yticks(np.arange(0, N, 10))

    # set minor ticks
    ax.set_xticks(np.arange(0, N, 5), minor=True)
    ax.set_yticks(np.arange(0, N, 5), minor=True)
    ax.grid(which="minor", alpha=0.15)
    ax.grid(alpha=0.4)

    # plot diagonal
    ax.plot([0, N], [0, N], color="black", lw=0.5, alpha=0.5, ls="--")


plt.tight_layout()
svfig("pairs_on_matrices_triu")
plt.show()
# %%

for M, lab in zip([Ms, Mb, Me], ["SNPs", "Blocks", "Edges"]):
    print(f"------------ {lab} ------------")
    # find positions of elements with value > 10
    wx, wy = np.where(M >= 14)
    for ix, iy in zip(wx, wy):
        if ix > iy:
            continue
        # print(ix, iy)
        s = f"{str_order[ix]:<20} {ix:>2} | {str_order[iy]:<20} {iy:>2} -- {M[ix, iy]:.0f}"
        print(s)

# %%
