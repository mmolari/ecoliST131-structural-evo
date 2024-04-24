# %%
import numpy as np
import pandas as pd
import pypangraph as pp
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import utils as ut

data_fld = pathlib.Path("data")
fig_fld = pathlib.Path("figs/f01")
fig_fld.mkdir(parents=True, exist_ok=True)

hdf = pd.read_csv(data_fld / "hotspots.csv", index_col=0)
hdf[["n_categories", "pangenome_len"]]

fname = "../../results/ST131_ABC/distances/summary-asm20-100-5.csv"
ddf = pd.read_csv(fname)
mask = ddf["si"] > ddf["sj"]
ddf = ddf[mask]
ddf.set_index(["si", "sj"], inplace=True)

hdf
# %%
# hs = "VZTFXIZVXB_f__YOCIMVGHSL_f"
# hs = "CGQZOLYSTD_r__ZLLQQUXUIP_r"
# hs = "JVNRLCFAVD_f__PLTCZQCVRD_r"
# hs = "IHKFSQQUKE_r__KPBYGJHRZJ_f"
hs = "GPKQYOCEJI_r__NKVSUZGURN_f"  # well-behaved
# hs = "XXVMWZCEKI_r__YUOECYBHUS_r"  # well-behaved
# hs = "BWEZXGGFBK_r__MVMOFPVELT_r"  # shortest
fname = f"../../results/ST131_ABC/backbone_joints/asm20-100-5/joints_pangraph/{hs}.json"
pan = pp.Pangraph.load_json(fname)


# %%
df = ut.extract_hotspot_stats(pan)

df["core_div_filtered"] = ddf.loc[df.index, "core_div_filtered"]

# %%

df["div"] = df["core_div_filtered"]
# df["div"] = df["local_core_div"]
# df = df[df["div"] < 0.00005]
fig, axs = plt.subplots(3, 3, figsize=(12, 12), sharex=True)

for nax, k in enumerate(
    [
        "n. private blocks",
        "merge_n_private",
        "len. private blocks",
        "n. breakpoints",
        "local_core_div",
        "merge_n_edges",
        "merge_core_len",
        "merge_acc_len",
        "merge_n_shared",
    ]
):
    ax = axs.flatten()[nax]
    sns.scatterplot(
        data=df,
        x="div",
        y=k,
        alpha=0.1,
        marker=".",
        ax=ax,
    )
    sns.lineplot(
        data=df,
        x="div",
        y=k,
        alpha=0.5,
        color="k",
        estimator=np.mean,
        # errorbar="sd",
        ax=ax,
    )

    # group by x in bins of size 0.0001
    df["bin"] = pd.cut(df["div"], 11)
    # x coordinate
    xc = df.groupby("bin", observed=True)["div"].mean()
    yc = df.groupby("bin", observed=True)[k].mean()
    err = df.groupby("bin", observed=True)[k].std()
    ax.errorbar(xc, yc, yerr=err, fmt="o", color="C1", ls="--")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / f"{hs}.png")
plt.show()

# %%
