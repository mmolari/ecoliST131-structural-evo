# %%

import numpy as np
import pandas as pd
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as sps


svfld = pathlib.Path("figs/1b")
svfld.mkdir(exist_ok=True, parents=True)

# %%
df = pd.read_csv("data/tips_df.csv", index_col=0)

df["has_acc"] = df["n_acc"] > 0
df["has_dupl"] = df["n_dupl"] > 0
df["has_event"] = df["n_events"] > 0
df
# %%


for i, col in enumerate(["n_acc", "n_dupl", "n_events"]):
    fig, axs = plt.subplots(2, 1, figsize=(5, 6), sharex=True)
    ax = axs[0]
    bins_x = np.linspace(0, df["branch_len"].max(), 15)
    bins_y = np.arange(df[col].max() + 2) - 0.5
    sns.histplot(data=df, x="branch_len", y=col, ax=ax, bins=(bins_x, bins_y))
    ax.set_xlabel("branch length")
    ylab = {
        "n_events": "n. terminal branch events",
        "n_dupl": "n. unique duplicated blocks",
        "n_acc": "n. unique accessory blocks",
    }[col]
    ax.set_ylabel(ylab)

    ax = axs[1]
    sns.regplot(
        data=df,
        x="branch_len",
        y={
            "n_events": "has_event",
            "n_dupl": "has_dupl",
            "n_acc": "has_acc",
        }[col],
        ax=ax,
        logistic=True,
        x_bins=bins_x,
        marker="o",
        # scatter_kws={"alpha": 0.3},
    )
    ax.set_xlabel("branch length")
    ax.set_ylabel("n > 0")

    sns.despine()
    plt.tight_layout()
    plt.savefig(svfld / f"branch_len_vs_{col}_logit.png", dpi=300)
    plt.show()

# %%

fig, axs = plt.subplots(2, 1, figsize=(5, 6), sharex=True)

ax = axs[0]
sns.histplot(data=df, x="branch_len", y="n_terminal_dupl_islands", bins=15, ax=ax)
ax.set_xlabel("branch length")
ax.set_ylabel("n. unique duplication islands")

ax = axs[1]
sns.regplot(
    data=df,
    x="branch_len",
    y="n_terminal_dupl_islands",
    scatter=True,
    scatter_kws={"alpha": 0.4, "color": "C0"},
    ax=ax,
)
# annotate p-value and correlation coefficient
r, p = sps.pearsonr(df["branch_len"], df["n_terminal_dupl_islands"])
ax.text(
    0.05,
    0.95,
    f"r = {r:.2f}\np = {p:.2e}",
    transform=ax.transAxes,
    ha="left",
    va="top",
)
ax.set_xlabel("branch length")
ax.set_ylabel("n. unique duplication islands")

plt.tight_layout()
plt.savefig(svfld / "dupl_islands.png", dpi=300)
plt.show()
# %%

fig, ax = plt.subplots(1, 1, figsize=(4, 3.5))

bins = np.arange(df[["n_acc", "n_events"]].max().max() + 2) - 0.5
sns.histplot(data=df, x="n_acc", y="n_events", bins=(bins, bins), ax=ax)
ax.set_xlabel("n. unique accessory blocks")
ax.set_ylabel("n. unique events")

plt.tight_layout()
plt.savefig(svfld / "acc_vs_events.png", dpi=300)
plt.show()

# %%
fig, ax = plt.subplots(1, 1, figsize=(4, 3.5))
sns.histplot(x=df["branch_len"], y=df["L_acc"] + df["L_dupl"], bins=15, ax=ax)
ax.set_xlabel("branch length")
ax.set_ylabel("private genome length (bp)")

plt.tight_layout()
plt.savefig(svfld / "private_genome_length.png", dpi=300)
plt.show()

# %%
