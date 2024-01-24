# %%
import pandas as pd
import subprocess
import pathlib
from Bio import Phylo
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


dset = "ST131_ABC"


fig_fld = pathlib.Path("figs/f06")
fig_fld.mkdir(exist_ok=True, parents=True)


# load joint coordinates dictionary
dset = "ST131_ABC"
fld = pathlib.Path(f"../../results/{dset}")

# %%
junction = "ANHYQUDQAA_r__IMAMHFLCWS_f"
# junction = "ATFHQYFPNW_f__GUDQDOMJFJ_f"
# junction = "JVNRLCFAVD_f__PLTCZQCVRD_r"
# junction = "BNIGRCVPML_f__BRZJUAZYZF_f" # IS
junction = "QHIDYNXVPV_f__QXRGVLICMC_r"  # IS


pan_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{junction}.json"

exp_fld = fig_fld / junction

# pangraph export
command = f"pangraph export {pan_file} -o {exp_fld} -ell 0 -nd -ll 0"
subprocess.run(command, shell=True, check=True)

# %%

# load tree
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
tree.ladderize()

# draw tree
fig, ax = plt.subplots(figsize=(3, 10))
Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda x: "")
ax.set_axis_off()
ax.plot([0, 1e-5], [100, 100], color="k", linewidth=1)
ax.text(0.5e-5, 99, "1e-5", ha="center", va="bottom", fontsize=10)

fig.savefig(fig_fld / "tree.svg", bbox_inches="tight", pad_inches=0)
plt.show()

# %%


# load junction info dataframes
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
df_j = pd.read_csv(df_file, index_col=0)
df_j["delta_len"] = df_j["max_length"] - df_j["min_length"]

# %%

# joint plot
fig, axs = plt.subplots(
    2,
    2,
    figsize=(5, 5),
    sharex="col",
    sharey="row",
    gridspec_kw={"width_ratios": [1, 0.2], "height_ratios": [0.2, 1]},
)

ax = axs[1, 0]
sns.scatterplot(
    data=df_j,
    x="delta_len",
    y="n_categories",
    alpha=0.2,
    ax=ax,
    rasterized=True,
)
ax.set_xlabel(r"$\Delta$ length (bp)")
ax.set_ylabel("n. different paths")

sns.histplot(
    data=df_j,
    x="delta_len",
    alpha=0.2,
    ax=axs[0, 0],
    legend=False,
    log_scale=True,
    bins=25,
    element="step",
)

sns.histplot(
    data=df_j,
    y="n_categories",
    alpha=0.2,
    ax=axs[1, 1],
    legend=False,
    log_scale=True,
    bins=25,
    element="step",
)

axs[0, 1].set_visible(False)
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "joint_plot.svg")
plt.show()

# %%
mask = np.abs(df_j["delta_len"] - 1000) < 200
mask &= df_j["n_categories"] == 2
df_j[mask]

# %%
