# %%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pathlib

fig_fld = pathlib.Path("figs/f00")
fig_fld.mkdir(parents=True, exist_ok=True)

# %%


def load_df():
    fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
    df = pd.read_csv(fname, index_col=0)

    fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/edge_pangenome.csv"
    df2 = pd.read_csv(fname, index_col=0)

    df = pd.merge(df, df2, on="edge", validate="one_to_one")
    df["delta_L"] = df["max_length"] - df["min_length"]
    df["acc_L"] = df["mean_length"] * 22
    return df


df = load_df()


def add_ann_info(df):
    fnames = {
        "df": "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/defensefinder_real.csv",
        "gm": "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/genomad_real.csv",
        "if": "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/integronfinder_real.csv",
        "is": "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/ISEScan_real.csv",
    }

    for k, fname in fnames.items():
        df2 = pd.read_csv(fname, index_col=0)
        df2 = df2["junction"].value_counts()
        df[f"{k}"] = df.index.map(df2)
        df[f"{k}"] = df[f"{k}"].fillna(0)
        df[f"{k}"] = df[f"{k}"].astype(int)

    return df


df = add_ann_info(df)


# %%
norm = mpl.colors.Normalize(vmin=0, vmax=1)
mapp = mpl.cm.ScalarMappable(norm=norm, cmap="coolwarm")

fig, ax = plt.subplots(1, 1, figsize=(8, 5))
sns.scatterplot(
    data=df,
    x="delta_L",
    y="n_categories",
    alpha=0.3,
    hue="nonempty_freq",
    palette="coolwarm",
    hue_norm=norm,
    legend=False,
)
plt.colorbar(mapp, ax=ax, label="occupation frequency")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"$\Delta L$ (bp)")
plt.ylabel("n. paths")
plt.tight_layout()
plt.savefig(fig_fld / "dL_vs_ncat_hue_freq.png")
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(8, 5))
sns.scatterplot(
    data=df,
    x="pangenome_len",
    y="n_categories",
    alpha=0.3,
    hue="nonempty_freq",
    palette="coolwarm",
    hue_norm=norm,
    legend=False,
)
plt.colorbar(mapp, ax=ax, label="occupation frequency")
plt.yscale("log")
plt.xscale("log")
plt.xlabel("pangenome length (bp)")
plt.ylabel("n. paths")
plt.tight_layout()
plt.savefig(fig_fld / "pL_vs_ncat_hue_freq.png")
plt.show()


# %%

info = {
    "is": ("insertion sequences", "C0", 4),
    "if": ("integrons", "C1", 5),
    "gm": ("prophages", "C2", 6),
    "df": ("defense islands", "C3", 7),
}

fig, ax = plt.subplots(1, 1, figsize=(7, 5))

x = df["delta_L"]
y = df["n_categories"]
ax.scatter(x, y, label=None, c="k", marker=".", alpha=0.3, zorder=-2)

for k, (label, c, m) in info.items():
    mask = df[k] > 0
    x = df[mask]["delta_L"]
    y = df[mask]["n_categories"]
    ax.scatter(x, y, label=label, c=c, marker=m, alpha=0.5)
plt.legend(title="annotations")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"$\Delta L$")
plt.ylabel("n. paths")
plt.tight_layout()
plt.savefig(fig_fld / "dL_vs_ncat_hue_ann.png")
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(7, 5))

x = df["pangenome_len"]
y = df["n_categories"]
ax.scatter(x, y, label=None, c="k", marker=".", alpha=0.3, zorder=-2)

for k, (label, c, m) in info.items():
    mask = df[k] > 0
    x = df[mask]["pangenome_len"]
    y = df[mask]["n_categories"]
    ax.scatter(x, y, label=label, c=c, marker=m, alpha=0.5)
plt.legend(title="annotations")
plt.yscale("log")
plt.xscale("log")
plt.xlabel("pangenome length (bp)")
plt.ylabel("n. paths")
plt.tight_layout()
plt.savefig(fig_fld / "pL_vs_ncat_hue_ann.png")
plt.show()


fig, ax = plt.subplots(1, 1, figsize=(7, 5))

x = df["acc_L"]
y = df["n_categories"]
ax.scatter(x, y, label=None, c="k", marker=".", alpha=0.3, zorder=-2)

for k, (label, c, m) in info.items():
    mask = df[k] > 0
    x = df[mask]["acc_L"]
    y = df[mask]["n_categories"]
    ax.scatter(x, y, label=label, c=c, marker=m, alpha=0.5)
plt.legend(title="annotations")
plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"$\Sigma \, L$ (bp)")
plt.ylabel("n. paths")
plt.tight_layout()
plt.savefig(fig_fld / "aL_vs_ncat_hue_ann.png")
plt.show()


# %%

fig, axs = plt.subplots(1, 2, figsize=(12, 5))

ax = axs[0]
sns.scatterplot(
    data=df,
    x="pangenome_len",
    y="pangenome_n_blocks",
    alpha=0.4,
    edgecolors="k",
    facecolors="none",
    ax=ax,
)
ax.set_xlabel("pangenome length (bp)")
ax.set_ylabel("n. blocks")

ax = axs[1]
sns.scatterplot(
    data=df,
    x="pangenome_len",
    y="delta_L",
    alpha=0.4,
    edgecolors="k",
    facecolors="none",
    ax=ax,
)
ax.set_xlabel("pangenome length (bp)")
ax.set_ylabel(r"$\Delta L$ (bp)")

for ax in axs:
    ax.set_yscale("log")
    ax.set_xscale("log")
plt.tight_layout()
plt.savefig(fig_fld / "pl_vs_dl.png")
plt.show()


# %%
sns.scatterplot(df, x="acc_L", y="n_categories", alpha=0.3, hue="nonempty_freq")
plt.yscale("log")
plt.xscale("log")

# %%
