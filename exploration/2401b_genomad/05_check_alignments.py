# %%
import pandas as pd
import pathlib
from Bio import SeqIO, Phylo
import subprocess
import pypangraph as pp
import numpy as np
import matplotlib.pyplot as plt

dset = "ST131_ABC"


fig_fld = pathlib.Path("figs/f05")
fig_fld.mkdir(exist_ok=True, parents=True)


data_fld = pathlib.Path("data")
pp_to_j = pd.read_csv(data_fld / "phages_to_joints.csv", index_col=0)

dset = "ST131_ABC"

# load joint coordinates dictionary
fld = pathlib.Path(f"../../results/{dset}")


# load tree
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
strains = [x.name for x in tree.get_terminals()]


# isolates genome length
iso_L = {}
for iso in strains:
    fa_fname = fld / f"../../data/fa/{iso}.fa"
    with open(fa_fname) as f:
        iso_L[iso] = len(SeqIO.read(f, "fasta"))

data_fld = pathlib.Path("data/phage_aln")
data_fld.mkdir(exist_ok=True, parents=True)

# %%
junction = "ANHYQUDQAA_r__IMAMHFLCWS_f"
# junction = "ATFHQYFPNW_f__GUDQDOMJFJ_f"

mask = pp_to_j["junction"] == junction

records = {}
for _, row in pp_to_j[mask].iterrows():
    iso = row["iso"]
    pb, pe = row["pb"], row["pe"]
    fa_fname = fld / f"../../data/fa/{iso}.fa"
    with open(fa_fname) as f:
        seq = SeqIO.read(f, "fasta")
    phage_seq = seq.seq[pb:pe]
    if not row["strand"]:
        phage_seq = phage_seq.reverse_complement()
    print(f"{iso} {row['strand']} {pb} {pe} {len(phage_seq)}")
    if iso in records:
        raise ValueError(f"Duplicate iso {iso}")
    records[iso] = SeqIO.SeqRecord(
        phage_seq, id=iso, description=f"{row['strand']} {pb} {pe}"
    )

records = [records[iso] for iso in strains if iso in records]
# save as fasta file
with open(data_fld / f"{junction}.fa", "w") as f:
    SeqIO.write(records, f, "fasta")

# %%

out_file = data_fld / f"{junction}.aln.fa"
if not out_file.exists():
    # align using mafft
    command = f"mafft --auto --thread 8 {data_fld / junction}.fa > {out_file}"
    subprocess.run(command, shell=True, check=True)

# visualize with aliview
command = f"aliview {out_file}"
subprocess.run(command, shell=True, check=True)

# %%

# load junction pangraph
pn_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{junction}.json"
pan = pp.Pangraph.load_json(pn_file)
bdf = pan.to_blockstats_df()
# %%


def aln_to_nogap_M(aln):
    N = len(aln)
    L = len(aln[0])
    M = np.empty((N, L), dtype="S1")
    for i, s in enumerate(aln):
        M[i, :] = np.array(list(s))
    has_gap = np.any(M == b"-", axis=0)
    M = M[:, ~has_gap]
    return M


def A_to_dist(A):
    "From alignment matrix to SNP distance matrix"
    N, L = A.shape
    D = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            D[i, j] = np.sum(A[i, :] != A[j, :])
    D += D.T
    return D / L


def order_strains(occs, strains):
    S = [o[0] for o in occs]
    S_ord = [s for s in strains if s in S]
    order = [S.index(s) for s in S_ord]
    return S_ord, order


mask = bdf["count"] > 10
mask &= bdf["len"] > 5000
mask &= ~bdf["duplicated"]

print(bdf[mask])

S_ords = {}

for bid in bdf[mask].index:
    block = pan.blocks[bid]
    aln, occs = block.alignment.generate_alignments()
    A = aln_to_nogap_M(aln)
    N, L = A.shape
    print(f"{bid} {N=} {L=}")
    D = A_to_dist(A)
    S_ord, order = order_strains(occs, strains)
    D = D[np.ix_(order, order)]
    S_ords[bid] = S_ord
    # plt.imshow(D, vmax=0.05, vmin=0)
    # plt.colorbar()
    # plt.title(f"{bid} {N=} {L=}")
    # plt.show()

# %%
col_file = f"figs/f04/J_structure_{junction}_colors.csv"
colors = pd.read_csv(col_file, index_col=0)

S = S_ords["BKVHJLHQTJ"]

# plot tree with only S isolates
S_tree = Phylo.read(tree_file, "newick")
for node in S_tree.get_terminals():
    if node.name not in S:
        S_tree.prune(node)
S_tree.ladderize()


fig, axs = plt.subplots(3, 2, figsize=(8, 10), sharey="row")
ax = axs[0, 0]
Phylo.draw(
    S_tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda x: ""
)
ax.grid(alpha=0.2)
ax.set_xlabel("position (bp)")

ax = axs[0, 1]
y = 0
for s in S:
    y += 1
    x = 0
    path = pan.paths[s]
    bids = path.block_ids
    for b in bids:
        l = bdf.loc[b, "len"]
        c = colors.loc[b, "Color"]
        ax.plot([x, x + l], [y, y], lw=3, color=c)
        x += l
ax.grid(alpha=0.2)

Bs = ["SVGDCUWZRJ", "AXINIMLYML", "BKVHJLHQTJ", "IZGMKIJGLP"]
for n, bid in enumerate(Bs):
    ax_x = n % 2
    ax_y = n // 2 + 1
    ax = axs[ax_y, ax_x]
    block = pan.blocks[bid]
    aln, occs = block.alignment.generate_alignments()
    A = aln_to_nogap_M(aln)
    N, L = A.shape
    print(f"{bid} {N=} {L=}")
    D = A_to_dist(A)
    S_ord, order = order_strains(occs, S)
    D = D[np.ix_(order, order)]
    ax.set_aspect("auto")
    im = ax.matshow(D, vmax=0.04, vmin=0)
    ticks = np.arange(len(S_ord))
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_title(f"{bid} {N=} {L=}", color=colors.loc[bid, "Color"], size=12)

plt.tight_layout()
plt.colorbar(im, ax=axs[1:, :], label="SNP density", shrink=0.5)
plt.savefig(fig_fld / f"J_{junction}_block_distance.png", dpi=300)
plt.show()

# %%

import itertools


def pairwise_divergence(A):
    N, L = A.shape
    D = np.zeros(L)
    for n1, n2 in itertools.combinations(range(N), 2):
        D += A[n1, :] != A[n2, :]
    return D / (N * (N - 1) / 2)


Div = {}

for n, bid in enumerate(Bs):
    block = pan.blocks[bid]
    aln, occs = block.alignment.generate_alignments()
    A = aln_to_nogap_M(aln)
    N, L = A.shape
    print(f"{bid} {N=} {L=}")
    S_ord, order = order_strains(occs, S)
    Div[bid] = pairwise_divergence(A)

# %%


def rolling_average(div, w):
    x = np.arange(len(div))
    y = div
    kernel = np.ones(w) / w
    x_avg = np.convolve(x, kernel, "valid")
    y_avg = np.convolve(y, kernel, "valid")
    return x_avg, y_avg


fig, axs = plt.subplots(len(Bs), 1, figsize=(7, 7), sharey=True, sharex=True)

for n, bid in enumerate(Bs):
    ax = axs[n]
    block = pan.blocks[bid]
    aln, occs = block.alignment.generate_alignments()
    A = aln_to_nogap_M(aln)
    N, L = A.shape
    S_ord, order = order_strains(occs, S)
    div = Div[bid]
    c = colors.loc[bid, "Color"]
    avg_pos, avg_div = rolling_average(div, 100)
    ax.plot(avg_pos, avg_div, label=bid, color=c)
    ax.set_ylabel("SNP density")
    ax.legend()


axs[-1].set_xlabel("position (bp)")
plt.tight_layout()
plt.savefig(fig_fld / f"J_{junction}_block_divergence.png", dpi=300)
plt.show()

# %%
