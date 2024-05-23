# %%
import pypangraph as pp
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from collections import defaultdict
import pathlib
import path_utils as pu

fg_fld = pathlib.Path("figs/n03")
fg_fld.mkdir(exist_ok=True, parents=True)

ref = "NZ_CP096110.1"
thr_len = 500

pan = pp.Pangraph.load_json("data/pangraph.json")
bdf = pan.to_blockstats_df()

# %%


def aln_vec(aln, strand):
    A = np.array([list(a) for a in aln])

    mask = np.any(A == "-", axis=0)
    A = A[:, ~mask]
    V = A[0] != A[1]
    if not strand:
        V = V[::-1]
    return V


# ========== Core Alignment ============

mask = bdf["core"] & (bdf["len"] >= thr_len)
CB = bdf[mask].index

CL = bdf[mask]["len"].sum()
print(f"core genome length = {CL/1e6} Mbp")

path = pan.paths[ref]
x = 0
breakpoints = [x]
muts = []
for i, b in enumerate(path.block_ids):
    if b not in CB:
        continue
    block = pan.blocks[b]
    strand = path.block_strands[i]
    aln, occ = block.alignment.generate_alignments()

    V = aln_vec(aln, strand)
    M = np.argwhere(V).flatten() + x
    x += len(V)
    breakpoints.append(x)
    muts += list(M)
# %%
fig, ax = plt.subplots(1, 1, figsize=(8, 3))
bw = 25000
bins = np.arange(0, breakpoints[-1] + bw, bw)
ax.hist(muts, bins=bins, weights=np.ones_like(muts) / bw)
# for bp in breakpoints:
#     ax.axvline(bp, alpha=0.1, lw=1, c="k")
ax.set_xlabel("core genome alignment")
ax.set_ylabel("fraction mutated sites")
ax.set_xlim(0, breakpoints[-1])
sns.despine()
plt.tight_layout()
plt.savefig(fg_fld / "mut_dens.png")
plt.savefig(fg_fld / "mut_dens.pdf")
plt.show()

# %%
# ============= merge paths =============

sdf = bdf.copy()
mask = sdf["core"] & (sdf["len"] > 500)
sdf = sdf[mask]

paths = pu.pangraph_to_path_dict(pan)
keep_f = lambda bid: bid in sdf.index
paths = pu.filter_paths(paths, keep_f)
Mg = pu.find_mergers(paths)

paths = pu.perform_mergers(Mg, paths, sdf)

paths["K12"].nodes = np.roll(paths["K12"].nodes, -2)

# %%

cmap = mpl.colormaps.get_cmap("nipy_spectral")
N = len(sdf)
block_colors = {b.id: cmap(i / (N - 1)) for i, b in enumerate(paths[ref].nodes)}

short_blocks = 2e4
fig, axs = plt.subplots(2, 1, figsize=(8, 5), sharex=True)

# length
ax = axs[0]
x = 0
for b in paths[ref].nodes:
    l = sdf.loc[b.id, "len"]
    c = block_colors[b.id]
    ax.bar(x, l, 1, color=c, edgecolor="k")
    x += 1
ax.set_yscale("log")
ax.set_ylabel("block length (bp)")
ax.axhline(short_blocks, zorder=-1, ls="--", color="silver")

# order
std_ord = {b.id: b.strand for b in paths[ref].nodes}
std_pos = {b.id: i for i, b in enumerate(paths[ref].nodes)}
ax = axs[1]
for lab, y in [("K12", 0), (ref, 1)]:
    path = paths[lab]
    x = 0
    for b in path.nodes:
        l = sdf.loc[b.id, "len"]
        c = block_colors[b.id]
        w = 0.6
        w = 0.4 if l < short_blocks else 0.6
        inv = b.strand != std_ord[b.id]
        ax.barh(y, 1, w, x, color=c, edgecolor=None)
        if inv:
            ax.arrow(
                x + 1,
                y - w / 2 - 0.2,
                -1,
                0,
                head_width=0.1,
                length_includes_head=True,
            )
        x += 1

ax.set_yticks([0, 1])
ax.set_yticklabels(["K12", ref])
ax.set_xlabel("minimal synteny blocks")

for x, b in enumerate(paths["K12"].nodes):
    c = block_colors[b.id]
    l = sdf.loc[b.id, "len"]
    if l < short_blocks:
        continue
    ax.plot([std_pos[b.id] + 0.5, x + 0.5], [0.7, 0.3], color=c)

sns.despine()
plt.tight_layout()
plt.savefig(fg_fld / "block_order.png", dpi=200)
plt.savefig(fg_fld / "block_order.svg")
plt.show()


# %%
def to_positions(P, B, S, O):
    N = len(B)
    assert N == len(P) - 1
    assert N == len(S)
    pos = defaultdict(list)
    for i in range(N):
        k = (P[i], P[(i + 1) % N], S[i], O[i])
        pos[B[i]].append(k)
    return pos


pos = {}
for strain in pan.strains():
    path = pan.paths[strain]
    P = path.block_positions
    O = path.block_nums
    B = path.block_ids
    S = path.block_strands
    pos[strain] = to_positions(P, B, S, O)

# %%

fig, ax = plt.subplots(1, 1, figsize=(6, 6))

mask = bdf["core"] & (bdf["len"] > thr_len)

for cb in bdf[mask].index:
    px = pos["K12"][cb]
    py = pos[ref][cb]
    assert len(px) == 1
    assert len(py) == 1
    px, py = px[0], py[0]

    x = px[:2]
    y = py[:2] if py[2] == px[2] else py[:2][::-1]
    if np.abs(np.diff(x)) > 1e6:
        continue
    color = block_colors[Mg[cb]] if cb in Mg else block_colors[cb]
    ax.plot(x, y, color=color)
ax.set_xlabel(f"K12 genome (bp)")
ax.set_ylabel(f"{ref} genome (bp)")
sns.despine()
plt.tight_layout()
plt.savefig(fg_fld / "block_order_dotplot.png", dpi=200)
plt.savefig(fg_fld / "block_order_dotplot.svg")
plt.show()

# %%
