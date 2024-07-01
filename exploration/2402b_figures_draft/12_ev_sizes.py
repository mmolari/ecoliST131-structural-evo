# %%

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib as pl

svfld = pl.Path("figs/f12")
svfld.mkdir(exist_ok=True, parents=True)

fname = "../../results/ST131_ABC/rates/asm20-100-5/terminal_coldspot.csv"
df = pd.read_csv(fname, index_col=0)

fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
jdf = pd.read_csv(fname, index_col=0)

dl = jdf["max_length"] - jdf["min_length"]
df["delta_l"] = dl
# %%

mask = df["event_type"] != "other"
df = df[mask]

fig, ax = plt.subplots()
sns.stripplot(
    data=df,
    x="event_type",
    y="delta_l",
    jitter=0.2,
    # color="black",
    # size=5,
    # alpha=0.5,
    ax=ax,
)
sns.boxplot(
    data=df,
    x="event_type",
    y="delta_l",
    showcaps=False,
    whis=[0, 100],
    width=0.6,
    color="silver",
    fliersize=0,
    linecolor="k",
    ax=ax,
)
plt.yscale("log")
plt.ylabel("Length difference")
plt.xlabel("Event type")
plt.tight_layout()
sns.despine()
plt.savefig(svfld / "delta_l.png")
plt.savefig(svfld / "delta_l.svg")
plt.show()


# %%

#  check alignments
import pypangraph as pp
from collections import Counter
import numpy as np

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
    D = G == 1
    I = G == N - 1

    res.append(
        {
            "block_id": bid,
            "length": L - 400,
            "length_all": L,
            "deletions_all": D.sum(),
            "insertions_all": I.sum(),
            "deletions": D[200:-200].sum(),  # remove 200 bp from each end
            "insertions": I[200:-200].sum(),
        }
    )
res = pd.DataFrame(res)
# %%

dcs = res["deletions"].cumsum()
ics = res["insertions"].cumsum()
lengths = res["length"].cumsum()

fig, ax = plt.subplots()
plt.plot(lengths, dcs, label="singleton deletion (bp)")
plt.plot(lengths, ics, label="singleton insertion (bp)")
plt.legend()
sns.despine()
plt.xlabel("core genome alignment")
plt.ylabel("cumulative sum")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(svfld / "singleton_indels.png")
plt.savefig(svfld / "singleton_indels.svg")
plt.show()
# %%

# same but for all positions
dcs = res["deletions_all"].cumsum()
ics = res["insertions_all"].cumsum()
lengths = res["length_all"].cumsum()

fig, ax = plt.subplots()
plt.plot(lengths, dcs, label="singleton deletion (bp)")
plt.plot(lengths, ics, label="singleton insertion (bp)")
plt.legend()
sns.despine()
plt.xlabel("core genome alignment")
plt.ylabel("cumulative sum")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(svfld / "singleton_indels_all.png")
plt.savefig(svfld / "singleton_indels_all.svg")
plt.show()
# %%
