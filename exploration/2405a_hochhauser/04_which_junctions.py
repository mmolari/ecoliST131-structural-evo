# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pathlib

fig_fld = pathlib.Path("figs/n04")
fig_fld.mkdir(parents=True, exist_ok=True)

odf = pd.read_csv("res/hs_junct_loc.csv")
mask = odf["in_acc"] != "no"
df = odf[mask]
df

# %%
jdf = pd.read_csv(
    "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
)
jdf_l = pd.read_csv(
    "../../results/ST131_ABC/backbone_joints/asm20-100-5/edge_pangenome.csv"
)
jdf = jdf.merge(jdf_l, on="edge")
jdf["in_hochhauser"] = jdf["edge"].isin(df["junction"])
jdf.set_index("edge", inplace=True)

# %%
freq_cmap = mpl.colormaps.get_cmap("coolwarm")
freq_norm = mpl.colors.Normalize(vmin=0, vmax=1)
len_cmap = mpl.colormaps.get_cmap("Greys")
len_norm = mpl.colors.LogNorm(vmin=1e2, vmax=1e5)
fig, ax = plt.subplots(1, 1, figsize=(2.5, 6))
for hs, gr in odf.groupby("id"):
    y = hs
    # print(gr)
    # ax.barh(y, len(gr), color="C0")
    mask = gr["in_acc"] != "no"
    # ax.barh(y, len(gr[mask]), color="C1")
    k = 0
    for i, row in gr[mask].iterrows():
        j = row["junction"]
        freq = jdf.loc[j, "nonempty_freq"]
        L = jdf.loc[j, "pangenome_len"]
        freq_color = freq_cmap(freq_norm(freq))
        len_color = len_cmap(len_norm(L))
        ax.scatter(k, y, color=freq_color, s=40, marker=4)
        ax.scatter(k, y, color=len_color, s=40, marker=5)
        k += 1
    if k == 0:
        ax.scatter(k, y, color="gray", s=40, marker="x")

ax.set_xticks([])
ax.set_ylabel("hotspot n.")

ax.set_xlim(-0.5, 6)
for k in ["top", "right", "bottom"]:
    ax.spines[k].set_visible(False)

plt.tight_layout()


# colorbars
cax1 = fig.add_axes([0.65, 0.15, 0.02, 0.3])
cax2 = fig.add_axes([0.65, 0.55, 0.02, 0.3])
cbar1 = mpl.colorbar.ColorbarBase(cax1, cmap=freq_cmap, norm=freq_norm)
cbar2 = mpl.colorbar.ColorbarBase(cax2, cmap=len_cmap, norm=len_norm)
cbar1.set_label("occupation frequency")
cbar2.set_label("pangenome length (bp)")

plt.savefig(fig_fld / "hs_junction_stats.png", dpi=200)
plt.show()

# %%
fig, ax = plt.subplots()
sns.scatterplot(
    data=jdf, x="pangenome_len", y="n_categories", alpha=0.5, ax=ax, hue="in_hochhauser"
)
ax.set_xscale("log")
ax.set_yscale("log")
plt.tight_layout()
plt.show()

# %%

np.random.seed(1)
fig, ax = plt.subplots()
sns.scatterplot(
    data=jdf, x="pangenome_len", y="n_categories", alpha=0.3, ax=ax, color="C0"
)
ylow = 0.025
dy = {}
for i, row in df.iterrows():
    j = row["junction"]
    hs_id = row["id"]
    x = jdf.loc[j, "pangenome_len"]
    y = jdf.loc[j, "n_categories"]
    ax.scatter(x, y, marker="o", color="red", facecolor="none")
    # random angle
    theta = np.random.uniform(0, 2 * np.pi)
    r = 8
    if hs_id not in dy:
        dy[hs_id] = ylow
        ylow += 0.03
    dx = 1.05
    ax.annotate(
        hs_id,
        (x, y),
        textcoords="axes fraction",
        xytext=(dx, dy[hs_id]),
        ha="center",
        va="center",
        fontsize=8,
        arrowprops=dict(
            arrowstyle="->",
            color="gray",
            lw=0.5,
            alpha=0.5,
            connectionstyle="angle,angleA=0,angleB=90,rad=15",
        ),
    )

ax.text(
    x=dx,
    y=ylow + 0.01,
    s="#HS",
    ha="center",
    va="center",
    fontsize=9,
    transform=ax.transAxes,
)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("pangenome length (bp)")
ax.set_ylabel("number of categories")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "hs_in_junctions_annotated.png", dpi=300)
plt.savefig(fig_fld / "hs_in_junctions_annotated.svg")
plt.show()

# %%
