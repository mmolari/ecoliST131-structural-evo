# %%
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import utils as ut
import pathlib

fig_fld = pathlib.Path("figs/n1")
fig_fld.mkdir(exist_ok=True, parents=True)

fname = "../../results/ST131_ABC/rates/asm20-100-5/terminal_coldspot.csv"
cdf = pd.read_csv(fname, index_col=0)
fname = "../../results/ST131_ABC/rates/asm20-100-5/terminal_branches.csv"
bdf = pd.read_csv(fname, index_col=0)
bdf["n_events"] = bdf["gain"] + bdf["loss"] + bdf["other"]


cdf = ut.assign_category(cdf)

# %%
colors = {
    "IS": "C0",
    "prophage": "C4",
    "integron": "C1",
    "none": "#b8b8b8",
}

fig, ax = plt.subplots(1, 1, figsize=(4, 5))
vc = cdf[["event_type", "cat"]].value_counts()

xlab = ["gain", "loss", "other"]
for x, k in enumerate(xlab):
    y = 0
    for c in ["IS", "prophage", "integron", "none"]:
        dy = vc.get((k, c), 0)
        kwargs = {"bottom": y, "color": colors[c]}
        if x == 0:
            kwargs["label"] = c
        ax.bar(x, dy, **kwargs)
        y += dy
    ax.text(x, y, f"{y}", ha="center", va="bottom")
ax.legend()

ax.set_xticks(range(len(xlab)))
ax.set_xticklabels(xlab)
ax.set_title("terminal coldspots")
ax.set_ylabel("n. events")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / f"terminal_coldspots.png")
plt.show()

# %%

fig, axs = plt.subplots(1, 2, figsize=(8, 4))

ax = axs[0]
sns.histplot(
    data=bdf,
    x="branch_length",
    y="n_events",
    discrete=(False, True),
    ax=ax,
)
ax.set_xlabel("branch length")
ax.set_ylabel("n. events")

ax = axs[1]
g = sns.histplot(
    data=bdf,
    x="branch_length",
    hue=bdf["n_events"] > 0,
    element="step",
    ax=ax,
)
ax.set_xlabel("branch length")
ax.set_ylabel("n. terminal branches")
# add legend
g.legend_.set_title("")
new_labels = ["no events", ">0 events"]
for t, l in zip(g.legend_.texts, new_labels):
    t.set_text(l)
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / f"terminal_branches.png")
plt.show()

# %%
