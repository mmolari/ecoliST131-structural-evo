# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# %%
df = pd.read_csv(
    "../../results/ST131/backbone_joints/asm20-100-5/junctions_stats.csv",
    index_col=0,
)
df["delta_len"] = df["max_length"] - df["min_length"]
df

# %%
fig, axs = plt.subplots(
    2,
    2,
    figsize=(8, 7),
    gridspec_kw={
        "height_ratios": [1.0, 0.3],
        "width_ratios": [0.3, 1.0],
    },
)

# scatterplot
ax = axs[0, 1]
sns.scatterplot(
    data=df,
    x="delta_len",
    y="n_categories",
    hue="has_dupl",
    ax=ax,
    alpha=0.3,
)
ax.set_xscale("symlog", linthresh=100)
ax.set_xlim(left=-50)
# ax.set_ylim(bottom=-0.5)
# ax.set_xlabel("max - min length (bp)")
# ax.set_ylabel("number of path categories")
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_yscale("log")
ax.legend(loc="upper left", title="has duplication")
ax.grid(alpha=0.3)

# marginal distribution length
ax = axs[1, 1]
sns.histplot(
    data=df,
    x="delta_len",
    ax=ax,
    element="step",
    stat="count",
    # hue="has_dupl",
    bins=[-10, 10] + list(np.logspace(1, 5, 20)),
)
ax.set_xscale("symlog", linthresh=100)
ax.set_xlim(left=-10)
ax.set_xlabel("max - min length (bp)")
ax.set_ylabel("count")
ax.grid(alpha=0.3)


# marginal distribution n_categories
ax = axs[0, 0]
sns.histplot(
    data=df,
    y="n_categories",
    ax=ax,
    element="step",
    stat="count",
    bins=np.arange(1, 100) - 0.5,
)
ax.set_yscale("log")
ax.set_xlabel("count")
ax.set_ylabel("number of path categories")
ax.grid(alpha=0.3)

# remove corner axis
ax = axs[1, 0]
ax.axis("off")


xlims = (-30, 1e6)
ylims = (0.8, 1e2)
axs[0, 1].set_xlim(*xlims)
axs[0, 1].set_ylim(*ylims)
axs[1, 1].set_xlim(*xlims)
axs[0, 0].set_ylim(*ylims)

sns.despine()
plt.tight_layout()
plt.savefig("figs/src_7/joints_overview.png", facecolor="white", dpi=300)
plt.savefig("figs/src_7/joints_overview.pdf")
plt.show()
# %%
