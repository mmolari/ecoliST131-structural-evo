# %%

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

import utils as ut

# %%

df = pd.read_csv(
    ut.expl_fld / "filtered_paths" / "block_stats_context.csv", index_col=0
)
svfld = ut.fig_fld / "context"
svfld.mkdir(exist_ok=True)

# %%
fig, axs = plt.subplots(1, 2, figsize=(12, 4))

ax = axs[0]
sns.histplot(
    data=df,
    x="count",
    hue="contexts",
    palette="tab10",
    multiple="stack",
    element="step",
    bins=np.arange(df["count"].max() + 1) + 0.5,
    ax=ax,
)

ax.set_xlabel("block frequency")
ax.set_ylabel("n. blocks")

ax = axs[1]
sns.histplot(
    data=df,
    x="len",
    hue="contexts",
    palette="tab10",
    multiple="stack",
    log_scale=True,
    ax=ax,
)

ax.set_xlabel("block length (bp)")
ax.set_ylabel("n. blocks")

sns.despine(fig)
plt.tight_layout()
fig.savefig(svfld / "9_context_hist.png")
plt.show()

# %%

sns.scatterplot(
    data=df,
    x="dist_intra",
    y="dist_inter",
    hue="contexts",
    palette="tab10",
    alpha=0.5,
)
plt.plot([0, 5e-5], [0, 5e-5], color="gray", ls="--")
plt.xlabel("avg. intra-context-clade distance")
plt.ylabel("avg. inter-context-clade distance")
sns.despine()
plt.tight_layout()
plt.savefig(svfld / "9_context_dist.png")
plt.show()

# %%

plt.hist2d(
    x=np.minimum(df["with"], df["without"]) - 1 + df["contexts"],
    y=df["gain"] + df["loss"] + df["move"],
    bins=(np.arange(0, 40, 1), np.arange(0, 40, 1)),
    norm=mpl.colors.LogNorm(),
    cmap="viridis",
)
plt.plot([0, 40], [0, 40], color="gray", ls="--")
plt.xlabel("maximal number of events")
plt.ylabel("inferred number of events")
plt.colorbar(label="n. blocks")
plt.tight_layout()
plt.savefig(svfld / "9_context_n_events.png")
plt.show()

# %%

pan = ut.load_pangraph()

# %%
b = pan.blocks["WNAPVOJCVV"]
# %%
b.sequence
# %%
