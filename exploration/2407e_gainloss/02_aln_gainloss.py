# %%

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib as pl
import pypangraph as pp
import numpy as np
import itertools as itt

svfld = pl.Path("figs/f2")
svfld.mkdir(exist_ok=True, parents=True)


fname = "../../results/ST131_ABC/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(fname)
bdf = pan.to_blockstats_df()
# %%
res = []
iso = "NZ_CP096110.1"
Bs = pan.paths[iso].block_ids
for k, bid in enumerate(Bs):
    if k % 10 == 0:
        print(f"{k}/{len(Bs)}", end="\r")
    b = pan.blocks[bid]
    if not bdf.loc[bid, "core"]:
        continue
    aln = b.alignment

    L = len(aln.consensus)
    if L < 500:
        continue
    A, O = b.alignment.generate_alignments()
    A = np.array([list(a) for a in A])
    N = A.shape[0]

    G = (A == "-").sum(axis=0)
    ter_deletions = G == 1
    ter_insertions = G == N - 1

    res.append(
        {
            "block_id": bid,
            "length": L,
            "ter_del": ter_deletions,
            "ter_ins": ter_insertions,
        }
    )
# %%
ter_modif = []
for r in res:
    bid = r["block_id"]
    L = r["length"]
    D = r["ter_del"]
    I = r["ter_ins"]

    for X, lab in zip([D, I], ["del", "ins"]):
        pos = 0
        for k, v in itt.groupby(X):
            l = len(list(v))
            if k:
                ter_modif.append(
                    {
                        "block_id": bid,
                        "block_length": L,
                        "length": l,
                        "pos": pos,
                        "type": lab,
                    }
                )
            pos += l
tdf = pd.DataFrame(ter_modif)

# %%
fig, axs = plt.subplots(
    2, 1, figsize=(6, 5), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
)
ax = axs[0]
sns.histplot(
    tdf,
    x="length",
    hue="type",
    bins=20,
    discrete=True,
    cumulative=True,
    element="step",
    fill=False,
    ax=ax,
)
ax.set_ylim(0, None)
ax.set_ylabel("n. in/dels")
ax.grid(axis="both", alpha=0.3)
ax = axs[1]
sns.stripplot(tdf, y="type", x="length", jitter=0.3, alpha=0.3, ax=ax)
ax.set_ylabel("")
ax.set_xlabel("in/del size (bp)")
plt.xscale("log")
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(svfld / "indel_size.png", dpi=300)
plt.show()

# %%

fig, ax = plt.subplots(figsize=(6, 3))
sns.histplot(
    tdf,
    x="length",
    hue="type",
    weights="length",
    cumulative=True,
    bins=20,
    discrete=True,
    element="step",
    fill=False,
    ax=ax,
)
ax.set_ylim(0, None)
plt.ylabel("cumulative in/del length (bp)")
plt.xlabel("in/del size (bp)")
plt.grid(axis="both", alpha=0.3)
plt.tight_layout()
plt.savefig(svfld / "indel_size_cumulative.png", dpi=300)
plt.show()


# %%
mask = (tdf["pos"] > 100) & (tdf["pos"] < tdf["block_length"] - 100)
fig, ax = plt.subplots(figsize=(6, 3))
sns.histplot(
    tdf[mask],
    x="length",
    hue="type",
    weights="length",
    cumulative=True,
    bins=20,
    discrete=True,
    element="step",
    fill=False,
    ax=ax,
)
ax.set_ylim(0, None)
plt.grid(axis="both", alpha=0.3)
plt.ylabel("cumulative in/del length (bp)")
plt.xlabel("in/del size (bp)")
plt.tight_layout()
plt.savefig(svfld / "indel_size_cumulative_midblock.png", dpi=300)
plt.show()

# %%
