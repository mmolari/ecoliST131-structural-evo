# %%
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import utils as ut
import pathlib

fig_fld = pathlib.Path("figs/n2")
fig_fld.mkdir(exist_ok=True, parents=True)

fname = "../../results/ST131_ABC/rates/asm20-100-5/nonsingleton_junct_df.csv"
cdf = pd.read_csv(fname, index_col=0)
fname = "../../results/ST131_ABC/rates/asm20-100-5/nonsingleton_branch_df.csv"
bdf = pd.read_csv(fname, index_col=0)
cdf["n_events"] = cdf["gain"] + cdf["loss"]

mask = ~cdf["undetermined"]
print(f"n. undetermined: {(~mask).sum()} / {mask.size}")
cdf = cdf[mask]

cdf = ut.assign_category(cdf)

# %%

fig, axs = plt.subplots(1, 2, figsize=(8, 4))

ax = axs[0]
sns.histplot(data=cdf, x="n_minority", y="n_events", ax=ax, discrete=True)
# plot diagonal
x = np.arange(0, cdf["n_minority"].max())
ax.plot(x, x, "k--", lw=1, c="gray")
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
# ax.set_aspect("equal")
ax.set_xlabel("n. minority")
ax.set_ylabel("n. events")

ax = axs[1]
colors = {
    "IS": "C0",
    "prophage": "C4",
    "integron": "C1",
    "none": "#b8b8b8",
}
for x, et in enumerate(["gain", "loss"]):
    y = 0
    for cat in cdf["cat"].unique():
        mask = cdf["cat"] == cat
        dy = cdf[mask][et].sum()
        kwargs = {
            "color": colors[cat],
            "width": 0.8,
        }
        if x == 0:
            kwargs["label"] = cat
        ax.bar(x, dy, bottom=y, **kwargs)
        y += dy
ax.legend()
ax.set_ylabel("n. events")
ax.set_xticks([0, 1])
ax.set_xticklabels(["gain", "loss"])

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "internal_gainloss.png")
plt.show()

# %%

fig, axs = plt.subplots(1, 2, figsize=(8, 4))

bdf["n_gainloss"] = bdf["n_gain"] + bdf["n_loss"]
mask = ~bdf["terminal"]

ax = axs[0]
sns.histplot(
    bdf[mask],
    x="branch_length",
    y="n_gainloss",
    ax=ax,
    discrete=(False, True),
)
ax.set_xlabel("branch length")
ax.set_ylabel("n. gain/loss events")

ax = axs[1]
g = sns.histplot(
    data=bdf[mask],
    x="branch_length",
    hue=bdf["n_gainloss"] > 0,
    element="step",
    stat="probability",
    common_norm=False,
    ax=ax,
)
ax.set_xlabel("branch length")
ax.set_ylabel("n. internal branches")
# add legend
g.legend_.set_title("")
new_labels = ["no events", ">0 events"]
for t, l in zip(g.legend_.texts, new_labels):
    t.set_text(l)

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "internal_branches.png")
plt.show()

# %%
