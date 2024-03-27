# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import pathlib
import gainloss_utils as ut
import json

fig_fld = pathlib.Path("figs/f05")
fig_fld.mkdir(exist_ok=True, parents=True)


def svfig(svname):
    pass
    plt.savefig(fig_fld / f"{svname}.png", dpi=200)
    plt.savefig(fig_fld / f"{svname}.svg")


edf = pd.read_csv(ut.fnames["events"])
fname = "../../results/ST131_ABC/rates/asm20-100-5/nonsingleton_branch_df.csv"
bdf = pd.read_csv(fname, index_col=0)

L_tot = bdf["branch_length"].sum()
L_ter = bdf[bdf["terminal"]]["branch_length"].sum()
L_int = bdf[~bdf["terminal"]]["branch_length"].sum()

# excluded junctions
j_count = edf["junction"].value_counts()
excl_j = j_count[j_count > 10].index

# %%

with open(
    "../../results/ST131_ABC/pangraph/asm20-100-5-alignment/filtered_corealignment_info_size.json"
) as f:
    aln_info = json.load(f)
core_aln_size = aln_info["core aln size"]
red_aln_size = aln_info["polished aln size"]

rates = []

for et in ["gain", "loss", "other"]:

    base_mask = edf["type"] == et

    stats = []

    # terminal singleton
    mask = base_mask & (edf["singleton"])
    stats.append((mask, "singleton", "terminal"))

    # terminal all
    mask = base_mask & (edf["terminal"])
    stats.append((mask, "all_raw", "terminal"))

    # internal
    mask = base_mask & (~edf["terminal"])
    stats.append((mask, "raw", "internal"))

    # all
    mask = base_mask
    stats.append((mask, "raw", "tot"))

    # terminal all corrected
    mask = base_mask & (edf["terminal"]) & (~edf["junction"].isin(excl_j))
    stats.append((mask, "all_corr", "terminal"))

    # internal corrected
    mask = base_mask & (~edf["terminal"]) & (~edf["junction"].isin(excl_j))
    stats.append((mask, "corr", "internal"))

    # all corrected
    mask = base_mask & (~edf["junction"].isin(excl_j))
    stats.append((mask, "corr", "tot"))

    for m, lab, L_lab in stats:
        match L_lab:
            case "terminal":
                L = L_ter
            case "internal":
                L = L_int
            case "tot":
                L = L_tot
            case _:
                raise ValueError(f"unknown L_lab: {L_lab}")

        n = m.sum()
        rates.append(
            {
                "type": et,
                "stat": L_lab + " | " + lab,
                "tree_part": L_lab,
                "n_events": n,
                "freq": L / n,
                "nmut_red": L * red_aln_size / n,
                "nmut_core": L * core_aln_size / n,
            }
        )

rdf = pd.DataFrame(rates)
rdf

# %%

df = rdf.pivot(index="type", columns="stat", values="nmut_red")
norm = mpl.colors.LogNorm(vmin=10, vmax=1000)
sns.heatmap(df, annot=True, fmt=".0f", cmap="coolwarm_r", norm=norm)
plt.title("n. (restricted) core-genome alignment mutations per event")
plt.xlabel("")
plt.ylabel("")
plt.tight_layout()
# plt.savefig(fig_fld / "nmut_red_heatmap.png", dpi=150)
plt.show()
# %%

tree_color = {
    "terminal": "C0",
    "internal": "C1",
    "tot": "C2",
}

for stat, plot_ylab in [
    ("freq", "tree length per event"),
    ("nmut_red", "n. restricted core-alignment mutations per event"),
    # ("nmut_core", "n. core-genome alignment mutations per event"),
]:

    fig, ax = plt.subplots(1, 1, figsize=(3, 4))
    # xlabs = ["gain", "loss", "other"]
    xlabs = ["gain", "loss"]
    for x, et in enumerate(xlabs):
        mask = rdf["type"] == et
        dx = -0.2
        for row in rdf[mask].itertuples():
            mec = None
            if row.stat.startswith("tot"):
                continue

            if "corr" in row.stat:
                mec = "k"
            elif "raw" in row.stat:
                mec = "gray"

            marker = "o"
            if "singleton" in row.stat:
                marker = "*"
            kwargs = {
                "color": tree_color[row.tree_part],
                "markeredgecolor": mec,
                "marker": marker,
            }
            y = getattr(row, stat)
            ax.plot([x + dx], [y], **kwargs)
            dx += 0.1

    # make legend
    for k, v in tree_color.items():
        if k == "tot":
            continue
        ax.plot([], [], "o", color=v, label=k)
    for k, mec in [("raw", "gray"), ("corrected", "k")]:
        ax.plot([], [], "o", color="white", markeredgecolor=mec, label=k)
    ax.plot([], [], "*", color="C0", label="singleton")
    ax.legend()

    if stat.startswith("nmut"):
        ax.set_ylim(10, 1000)

    ax.set_xticks(range(len(xlabs)))
    ax.set_xticklabels(xlabs)
    ax.set_yscale("log")
    ax.set_ylabel(plot_ylab)
    ax.grid(axis="y", which="major", linewidth=0.5, alpha=1)
    ax.grid(axis="y", which="minor", linewidth=0.5, alpha=0.5)
    sns.despine()
    plt.tight_layout()
    svfig(f"rate_{stat}")
    plt.show()


# %%
