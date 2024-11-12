# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib

res_fld = pathlib.Path("res")
res_fld.mkdir(exist_ok=True)

fig_fld = pathlib.Path("figs")
fig_fld.mkdir(exist_ok=True)

# load gene counts
Mc = pd.read_csv(res_fld / "mt_counts.csv", index_col=[0, 1, 2])  # mutation counts
Tc = pd.read_csv(res_fld / "tot_counts.csv", index_col=0).T  # total counts
Bc = pd.read_csv(res_fld / "bck_counts.csv", index_col=[0, 1, 2])  # background counts

# %% mutations filter
# only select genes with up to 3 mutations per kbp
Tc["select"] = (Tc["npol"] / Tc["ncol"]) <= 3 / 1000

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
sns.histplot(
    data=Tc,
    bins=100,
    x="ncol",
    y="npol",
    ax=ax,
    discrete=(False, True),
    hue="select",
)
plt.tight_layout()
plt.savefig(fig_fld / "mut_filter.png", dpi=300)
plt.show()

print(Tc["select"].mean())
selected_genes = Tc[Tc["select"]].index
Mc = Mc[selected_genes]
Bc = Bc[selected_genes]

# %%


def to_count_df(df):
    tot = (
        df.sum(axis=1)
        .reset_index()
        .rename(
            columns={0: "counts", "level_0": "from", "level_1": "to", "level_2": "syn"}
        )
    )

    # synonim counts
    sc = tot[tot["syn"] == "S"].set_index(["from", "to"])["counts"]
    # non-synonim counts
    nc = tot[tot["syn"] == "N"].set_index(["from", "to"])["counts"]

    return sc, nc


msc, mnc = to_count_df(Mc)
bsc, bnc = to_count_df(Bc)

cdf = pd.DataFrame(
    {
        "Sm": msc,
        "Nm": mnc,
        "Sb": bsc,
        "Nb": bnc,
    }
)

cdf["S_obs"] = cdf["Sm"] / cdf["Sb"]
cdf["N_obs"] = cdf["Nm"] / cdf["Nb"]
cdf["mut"] = cdf.index.map(lambda x: f"{x[0]}->{x[1]}")
cdf
# %%


def heatmap(S, N):
    fig, axs = plt.subplots(1, 2, figsize=(8, 3.8))
    ax = axs[0]
    sns.heatmap(S, ax=ax, cmap="viridis")
    ax = axs[1]
    sns.heatmap(N, ax=ax, cmap="viridis")
    return fig, axs


S = cdf.pivot_table(columns="to", index="from", values="S_obs")
N = cdf.pivot_table(columns="to", index="from", values="N_obs")

# heatmap for both mutations
fig, axs = heatmap(S, N)
axs[0].set_title("synonymous mutations rate")
axs[1].set_title("non-synonymous mutations rate")
plt.tight_layout()
plt.savefig(fig_fld / "mut_heatmap.png", dpi=300)
plt.show()

# %%


def errorbar_size(N, k):
    p = k / N
    # z = 1.96
    # z2 = z**2
    # num = p + z2 / (2 * N)
    # den = 1 + z2 / N
    # factor = z / np.sqrt(N)
    # return factor * np.sqrt(num * (1 - num) / den)
    return np.sqrt(p * (1 - p) / N)


cS = "#1f77b4"
cN = "#ff7f6e"

fig, axs = plt.subplots(
    3,
    1,
    figsize=(8, 8),
    sharex=True,
    height_ratios=[2, 1, 1],
)
ax = axs[0]
for i, (idx, row) in enumerate(cdf.iterrows()):
    # plot ratio observed
    ax.plot([i, i], [row["S_obs"], row["N_obs"]], color="black", alpha=0.5)
    ax.errorbar(
        i,
        row["S_obs"],
        yerr=errorbar_size(row["Sb"], row["Sm"]),
        fmt="o",
        color=cS,
        capsize=5,
        label="synonymous" if i == 0 else None,
    )
    ax.errorbar(
        i,
        row["N_obs"],
        yerr=errorbar_size(row["Nb"], row["Nm"]),
        fmt="o",
        color=cN,
        capsize=5,
        label="non-synonymous" if i == 0 else None,
    )
ax.set_ylabel("mut. rate (observed/possible)")
ax.set_yscale("log")
ax.set_ylim(bottom=1e-6)
ax.grid(alpha=0.3)
ax.legend()

for ax, color, key, lab in zip(
    axs[1:], [cS, cN], ["S", "N"], ["synonymous", "non-synonymous"]
):
    ax.bar(
        np.arange(len(cdf)) - 0.2,
        cdf[f"{key}b"],
        width=0.4,
        color=color,
        label="n. possible",
    )
    ax.bar(
        np.arange(len(cdf)) + 0.2,
        cdf[f"{key}m"],
        width=0.4,
        color=color,
        alpha=0.5,
        label="n. observed",
    )
    ax.legend()
    ax.set_yscale("log")
    ax.grid(alpha=0.3)
    ax.set_ylabel(f"n. {lab} muts")


ax.set_xticks(np.arange(len(cdf)))
ax.set_xticklabels(cdf["mut"], rotation=90)


plt.tight_layout()
plt.savefig(fig_fld / "mut_ratios.png", dpi=300)
plt.show()
# %%


fig, ax = plt.subplots(1, 1, figsize=(8, 4))
for i, (idx, row) in enumerate(cdf.iterrows()):
    # plot ratio observed
    ax.plot([i, i], [row["S_obs"], row["N_obs"]], color="black", alpha=0.5)
    ax.errorbar(
        i,
        row["S_obs"],
        yerr=errorbar_size(row["Sb"], row["Sm"]),
        fmt="o",
        color=cS,
        capsize=5,
        label="synonymous" if i == 0 else None,
    )
    ax.errorbar(
        i,
        row["N_obs"],
        yerr=errorbar_size(row["Nb"], row["Nm"]),
        fmt="o",
        color=cN,
        capsize=5,
        label="non-synonymous" if i == 0 else None,
    )
ax.set_ylabel("mut. rate (observed/possible)")
ax.set_yscale("log")
ax.set_ylim(bottom=1e-6)
ax.grid(alpha=0.3)
ax.legend()
ax.set_xticks(np.arange(len(cdf)))
ax.set_xticklabels(cdf["mut"], rotation=90)


plt.tight_layout()
plt.savefig(fig_fld / "mut_ratios_simple.png", dpi=300)
plt.show()

# %%

# total ratio:
sum_df = cdf.sum()
dn = sum_df["Nm"] / sum_df["Nb"]
ds = sum_df["Sm"] / sum_df["Sb"]
dnds = dn / ds

print(f"{dn=:.2e}, {ds=:.2e}, {dnds=:.2e}")
# %%
