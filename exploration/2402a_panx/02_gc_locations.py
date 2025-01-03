# %%
import pandas as pd
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
import utils as ut
from scipy.stats import entropy

fig_fld = pathlib.Path("figs/f02")
fig_fld.mkdir(parents=True, exist_ok=True)

gene_cluster_json = "../../data/panX/data/ST131_ABC/vis/geneCluster.json"
gc = ut.load_genecluster_json(gene_cluster_json)
gids = gc.index.to_numpy()

loc_fld_real = pathlib.Path("../../results/ST131_ABC/panx/gc_j_pos/asm20-100-5/real")
loc_fld_rand = pathlib.Path("../../results/ST131_ABC/panx/gc_j_pos/asm20-100-5/rand")


def parse_locations(gids, path):
    gc_info = []

    for gid in gids:
        res = {}
        row = gc.loc[gid]
        res["gid"] = gid
        res["n_iso"] = row["count"]
        res["n_loci"] = len(row["locus"].split())
        res["dupl"] = row["dupli"] == "yes"
        loc_file = path / f"genecl_{gid}.csv"
        if not loc_file.exists():
            # print(f"missing {loc_file}")
            raise ValueError

        try:
            df = pd.read_csv(loc_file)
            jc = df["junction"].value_counts()
            res["max_j"] = jc.max()
            res["j_nuniq"] = len(jc)
            res["j_entr"] = entropy(jc)

            res["acc"] = len(df)
            res["core"] = res["n_loci"] - len(df)
        except pd.errors.EmptyDataError:
            # file is empty
            # print(f"empty {loc_file}")
            res["core"] = res["n_loci"]

        res["ann"] = row["ann"]
        gc_info.append(res)
    gc_info = pd.DataFrame(gc_info)
    gc_info.set_index("gid", inplace=True)
    gc_info.fillna({"acc": 0, "max_j": 0, "j_nuniq": 0}, inplace=True)
    gc_info = gc_info.astype({"acc": int, "max_j": int, "j_nuniq": int})
    return gc_info


# %%

print("parsing real positions")
gci = parse_locations(gids, loc_fld_real)

print("parsing random positions")
gci_rand = parse_locations(gids, loc_fld_rand)
# %%
M = max(gci_rand["j_nuniq"].max(), gci["j_nuniq"].max())
bins = np.arange(0, M + 1, 5)

sns.histplot(gci, x="j_nuniq", bins=bins, element="step", label="real")
sns.histplot(gci_rand, x="j_nuniq", bins=bins, element="step", label="random")
plt.legend()
plt.tight_layout()
plt.show()

# %%

# check number of perfectly placed core genes

core_mask = gci["dupl"] == False
core_mask &= gci["n_loci"] == 222
gci[core_mask]

fig, axs = plt.subplots(1, 2, figsize=(8, 4))
ax = axs[0]
sns.histplot(
    gci[core_mask]["j_nuniq"], discrete=True, ax=ax, element="step", color="C0"
)
ax.set_title("core genes - real annotations")

ax = axs[1]
sns.histplot(
    gci_rand[core_mask]["j_nuniq"],
    discrete=True,
    ax=ax,
    element="step",
    color="C1",
)
ax.set_title("core genes - random annotations")

for ax in axs:
    ax.set_xlabel("n. of associated junctions")
    ax.set_ylabel("n. of gene clusters")

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "core_gc_junctions.png")
plt.show()
# %%
weird_mask = core_mask & (gci["j_nuniq"] == 1)
sns.histplot(gci[weird_mask]["acc"], discrete=True, element="step")
wid = gci[weird_mask].index

for gid in wid:
    print(gc.loc[gid]["locus"])
    print(gc.loc[gid]["ann"])
    print()
    pf = f"../../results/ST131_ABC/panx/gc_j_pos/asm20-100-5/real/genecl_{gid}.csv"
    df = pd.read_csv(pf)
    fig, ax = plt.subplots(1, 1, figsize=(12, 20))
    for i, row in df.iterrows():
        As, Ae = row["jab"], row["jae"]
        Gb, Ge = row["ib"], row["ie"]

        s = row["j_strand"]
        if not s:
            S, E = 0, As - Ae
            Gb, Ge = -Ge + Ae, -Gb + Ae
        else:
            S, E = As - As, Ae - Ae
            Gb, Ge = Gb - As, Ge - As
        ax.plot([S, E], [i, i], "k", lw=4, alpha=0.5)
        ax.plot([Gb, Ge], [i, i], "-|", color="r")
    ax.set_title(f"Gene cluster {gid}")
    plt.show()
    break

# %%

intermediate_mask = gci["dupl"] == False
intermediate_mask &= gci["n_loci"] < 222
gci[intermediate_mask]

fig, axs = plt.subplots(1, 2, figsize=(8, 4))
ax = axs[0]
sns.histplot(
    gci[intermediate_mask]["j_nuniq"], discrete=True, ax=ax, element="step", color="C0"
)
ax.set_title("non-core genes - real annotations")

ax = axs[1]
sns.histplot(
    gci_rand[intermediate_mask]["j_nuniq"],
    discrete=True,
    ax=ax,
    element="step",
    color="C1",
)
ax.set_title("non-core genes - random annotations")

for ax in axs:
    ax.set_xlabel("n. of associated junctions")
    ax.set_ylabel("n. of gene clusters")

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "intermediate_gc_junctions.png")
plt.show()

# %%

# %%


dupl_mask = gci["dupl"] == True
gci[dupl_mask]

fig, axs = plt.subplots(1, 2, figsize=(8, 4))
ax = axs[0]
sns.histplot(
    gci[dupl_mask]["j_nuniq"], discrete=True, ax=ax, element="step", color="C0"
)
ax.set_title("core genes - real annotations")

ax = axs[1]
sns.histplot(
    gci_rand[dupl_mask]["j_nuniq"],
    discrete=True,
    ax=ax,
    element="step",
    color="C1",
)
ax.set_title("core genes - random annotations")

for ax in axs:
    ax.set_xlabel("n. of associated junctions")
    ax.set_ylabel("n. of gene clusters")

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "duplicated_gc_junctions.png")
plt.show()


# %%


gci[dupl_mask].sort_values("j_nuniq", ascending=False).head(10)["ann"].iloc[4]

# never present in a core region
# mostly transposases
# %%
