# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib

fig_fld = pathlib.Path("figs/f01")
fig_fld.mkdir(parents=True, exist_ok=True)


fname = "../../results/ST131_ABC/rates/asm20-100-5/merged_events_df.csv"
edf = pd.read_csv(fname, index_col=0)
edf.head()

fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
jdf = pd.read_csv(fname, index_col=0)

dL = jdf["max_length"] - jdf["min_length"]
edf["dL"] = dL

# %%
mask = edf["terminal"]
bins = np.logspace(2, np.log10(edf["dL"].max() + 2), 40)
g = sns.FacetGrid(
    edf[mask],
    row="type",
    row_order=["gain", "loss"],
    hue="mge_cat",
    hue_order=[
        "IS",
        "none",
        "defense",
        "integron",
        "prophage",
    ],
    height=2.1,
    aspect=3,
).map(
    sns.histplot,
    "dL",
    bins=bins,
)

# add small-size legend to first plot
g.axes[0, 0].legend(title="MGE category", loc="upper right", fontsize="small")

plt.yscale("log")
plt.xscale("log")
plt.tight_layout()
plt.savefig(fig_fld / "dL_terminal_gainloss.png")
plt.show()

# %%
mask = edf["terminal"] & (edf["mge_cat"] == "IS")
edf[mask].sort_values("dL", ascending=False).head()

# %%
