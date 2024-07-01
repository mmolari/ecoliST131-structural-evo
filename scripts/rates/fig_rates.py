# %%
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pathlib
import argparse
import merge_event_df as mgut


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--coldspots_df", type=str)
    parser.add_argument("--events_df", type=str)
    parser.add_argument("--terminal_j_df", type=str)
    parser.add_argument("--internal_j_df", type=str)
    parser.add_argument("--terminal_b_df", type=str)
    parser.add_argument("--internal_b_df", type=str)
    parser.add_argument("--alignment_info", type=str)
    parser.add_argument("--out_fld", type=str)
    return parser.parse_args()


def load_data(args):

    with open(args.alignment_info, "r") as f:
        aln_info = json.load(f)
    aln_L = aln_info["polished aln size"]
    aln_L_core = aln_info["core aln size"]

    js = pd.read_csv(args.coldspots_df, index_col=0)
    js = mgut.assign_mge_category(js)

    edf = pd.read_csv(args.events_df)

    j_dfs = {}
    for k in ["terminal", "internal"]:
        j_dfs[k] = pd.read_csv(getattr(args, f"{k}_j_df"), index_col=0)
        j_dfs[k] = mgut.assign_mge_category(j_dfs[k])
    b_dfs = {}

    for k in ["terminal", "internal"]:
        b_dfs[k] = pd.read_csv(getattr(args, f"{k}_b_df"), index_col=0)
    k = "terminal"
    b_dfs[k]["n_events"] = b_dfs[k]["gain"] + b_dfs[k]["loss"] + b_dfs[k]["other"]
    b_dfs[k]["n_muts"] = b_dfs[k]["branch_length"] * aln_L
    k = "internal"
    b_dfs[k]["n_events"] = b_dfs[k]["n_gain"] + b_dfs[k]["n_loss"] + b_dfs[k]["n_other"]
    b_dfs[k]["n_muts"] = b_dfs[k]["branch_length"] * aln_L

    return aln_L, aln_L_core, js, j_dfs, b_dfs, edf


cat_colors = {
    "IS": "C0",
    "prophage": "C4",
    "integron": "C1",
    "none": "#b8b8b8",
    "defense": "C2",
}


def plot_coldspots_breakdown(js, fig_fld):

    fig, axs = plt.subplots(
        1,
        2,
        figsize=(7, 3.5),
        gridspec_kw={"width_ratios": [1, 2.3]},
        sharey=True,
    )

    # singletons
    ax = axs[0]
    colors = {
        True: "gray",
        False: "lightgray",
    }

    et = js["singleton"].value_counts().to_dict()
    xticks, xlabels = [], []
    x = 0
    for k, v in et.items():
        b = ax.bar(x, v, color=colors[k])
        ax.bar_label(b)
        xticks.append(x)
        xlabels.append("singleton" if k else "non-singleton")
        x += 1
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    ax.set_ylabel("n. binary junctions")

    ax = axs[1]
    et = js["cat"].value_counts().sort_values(ascending=False)
    x = 0
    xticks, xlabels = [], []
    for k, v in et.items():
        b = ax.bar(x, v, color=cat_colors[k])
        ax.bar_label(b)
        xticks.append(x)
        xlabels.append(k)
        x += 1

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)
    # ax.set_ylabel("n. binary junctions")

    sns.despine()
    plt.tight_layout()
    plt.savefig(fig_fld / "suppl_backbone_coldspots_breakdonw.pdf")
    plt.close()


def fig_terminal_coldspots_counts(cdf, fig_fld):

    def barplot_events(cdf, ax):
        vc = cdf[["event_type", "cat"]].value_counts()

        xlab = ["gain", "loss", "other"]

        for x, k in enumerate(xlab):
            y = 0
            for c in cdf["cat"].unique().sort_values():
                dy = vc.get((k, c), 0)
                kwargs = {"bottom": y, "color": cat_colors[c]}
                if x == 0:
                    kwargs["label"] = c
                ax.bar(x, dy, **kwargs)
                y += dy
            ax.text(x, y, y, ha="center", va="bottom")
        ax.legend()

        ax.set_xticks(range(len(xlab)))
        ax.set_xticklabels(xlab)
        ax.set_ylabel("n. events")

    fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.5))
    barplot_events(cdf, ax)
    sns.despine()
    plt.tight_layout()
    plt.savefig(fig_fld / "terminal_coldspots_count.pdf")
    plt.close()


def fig_terminal_branches_events(bdf, fig_fld):
    fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

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
    ax.set_ylabel("n. branches")
    # add legend
    g.legend_.set_title("")
    new_labels = ["no events", r"$\geq 1$" + " events"]
    for t, l in zip(g.legend_.texts, new_labels):
        t.set_text(l)
    sns.despine()
    plt.tight_layout()
    plt.savefig(fig_fld / "suppl_terminal_branches_vs_events.pdf")
    plt.close()


def fig_supple_terminal_muts_vs_events(bdf, fig_fld):
    fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

    ax = axs[0]
    sns.histplot(
        data=bdf,
        x="n_muts",
        y="n_events",
        discrete=(False, True),
        ax=ax,
    )
    ax.set_xlabel("n. mutations")
    ax.set_ylabel("n. events")

    ax = axs[1]
    g = sns.histplot(
        data=bdf,
        x="n_muts",
        hue=bdf["n_events"] > 0,
        element="step",
        ax=ax,
    )
    ax.set_xlabel("n. mutations")
    ax.set_ylabel("n. branches")
    # add legend
    g.legend_.set_title("")
    new_labels = ["no events", r"$\geq 1$" + " events"]
    for t, l in zip(g.legend_.texts, new_labels):
        t.set_text(l)
    sns.despine()
    plt.tight_layout()
    plt.savefig(fig_fld / "suppl_terminal_muts_vs_events.pdf")
    plt.close()


def fig_internal_events(cdf, fig_fld):

    fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

    ax = axs[0]
    g = sns.histplot(data=cdf, x="n_minority", y="n_events", ax=ax, discrete=True)

    # plot diagonal
    # x = np.arange(0, cdf["n_minority"].max())
    x = np.arange(0, cdf["n_events"].max())
    ax.set_yticks(range(0, x.max() + 1 + 2, 2))
    ax.plot(x, x, "--", lw=1, c="gray", label="diagonal")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    # ax.set_aspect("equal")
    ax.set_xlabel("n. isolates with minority pattern")
    ax.set_ylabel("n. events")
    # ax.legend(loc="center right")

    ax = axs[1]
    xlabs = ["gain", "loss", "other"]
    for x, et in enumerate(xlabs):
        y = 0
        for cat in cdf["cat"].unique().sort_values():
            mask = cdf["cat"] == cat
            dy = cdf[mask][et].sum()
            kwargs = {
                "color": cat_colors[cat],
                "width": 0.8,
            }
            if x == 0:
                kwargs["label"] = cat
            ax.bar(x, dy, bottom=y, **kwargs)
            y += dy
        ax.text(x, y, y, ha="center", va="bottom")
    ax.legend()
    ax.set_ylabel("n. events")
    ax.set_xticks(range(len(xlabs)))
    ax.set_xticklabels(xlabs)

    sns.despine()
    plt.tight_layout()

    # colorbar in inset
    cax = fig.add_axes([0.40, 0.45, 0.02, 0.35])
    plt.colorbar(
        g.collections[0],
        cax=cax,
        orientation="vertical",
        label="n. non-singleton junctions",
    )

    plt.savefig(fig_fld / "suppl_internal_gainloss.pdf")
    plt.close()


def fig_internal_branches(df, xlab, xkey, fname, fig_fld):
    fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

    mask = ~df["terminal"]
    sdf = df[mask]

    ax = axs[0]
    sns.histplot(
        sdf,
        x=xkey,
        y="n_events",
        ax=ax,
        discrete=(False, True),
    )
    ax.set_xlabel(xlab)
    ax.set_ylabel("n. events")

    ax = axs[1]
    g = sns.histplot(
        data=sdf,
        x=xkey,
        hue=sdf["n_events"] > 0,
        element="step",
        ax=ax,
    )
    ax.set_xlabel(xlab)
    ax.set_ylabel("n. branches")
    # add legend
    g.legend_.set_title("")
    new_labels = ["no events", r"$\geq 1$" + " events"]
    for t, l in zip(g.legend_.texts, new_labels):
        t.set_text(l)

    sns.despine()
    plt.tight_layout()
    plt.savefig(fig_fld / fname)
    plt.close()


def evaluate_rates(edf, bdf, red_aln_size, core_aln_size):

    L_tot = bdf["branch_length"].sum()
    L_ter = bdf[bdf["terminal"]]["branch_length"].sum()
    L_int = bdf[~bdf["terminal"]]["branch_length"].sum()

    j_count = edf["junction"].value_counts()
    excl_j = j_count[j_count > 10].index

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

        L_dict = {
            "terminal": L_ter,
            "internal": L_int,
            "tot": L_tot,
        }

        for m, lab, L_lab in stats:
            L = L_dict[L_lab]

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
                    "n_ev_per_mut_red": n / (L * red_aln_size),
                    "n_ev_per_mut_core": n / (L * core_aln_size),
                }
            )

    rdf = pd.DataFrame(rates)
    return rdf


def fig_rates_overview(rates, fig_fld):
    df = rates.pivot(index="type", columns="stat", values="nmut_red")
    norm = mpl.colors.LogNorm(vmin=10, vmax=1000)
    sns.heatmap(df, annot=True, fmt=".0f", cmap="coolwarm_r", norm=norm)
    plt.title("n. (restricted) core-genome alignment mutations per event")
    plt.xlabel("")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(fig_fld / "rates_heatmap.pdf")
    plt.close()


def fig_rates(rates, fig_fld):

    tree_color = {
        "terminal": "C0",
        "internal": "C1",
        "tot": "C2",
    }

    for stat, plot_ylab in [
        ("freq", "tree length per event"),
        ("nmut_red", "n. restricted core-alignment mutations per event"),
        ("n_ev_per_mut_red", "n. events per core-alignment mutation"),
    ]:

        fig, ax = plt.subplots(1, 1, figsize=(3, 4))
        # xlabs = ["gain", "loss", "other"]
        xlabs = ["gain", "loss"]
        for x, et in enumerate(xlabs):
            mask = rates["type"] == et
            dx = -0.2
            for row in rates[mask].itertuples():
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
        elif stat.startswith("n_ev"):
            ax.set_ylim(0.001, 0.1)

        ax.set_xticks(range(len(xlabs)))
        ax.set_xticklabels(xlabs)
        ax.set_yscale("log")
        ax.set_ylabel(plot_ylab)
        ax.grid(axis="y", which="major", linewidth=0.5, alpha=1)
        ax.grid(axis="y", which="minor", linewidth=0.5, alpha=0.5)
        sns.despine()
        plt.tight_layout()
        plt.savefig(fig_fld / f"rate_{stat}.pdf")
        plt.close()


if __name__ == "__main__":

    args = parse_args()

    fig_fld = pathlib.Path(args.out_fld)
    fig_fld.mkdir(exist_ok=True, parents=True)

    aln_L, aln_L_core, js, j_dfs, b_dfs, edf = load_data(args)

    plot_coldspots_breakdown(js, fig_fld)
    fig_terminal_coldspots_counts(j_dfs["terminal"], fig_fld)
    fig_terminal_branches_events(b_dfs["terminal"], fig_fld)
    fig_supple_terminal_muts_vs_events(b_dfs["terminal"], fig_fld)
    fig_internal_events(j_dfs["internal"], fig_fld)
    fig_internal_branches(
        b_dfs["internal"],
        "branch length",
        "branch_length",
        "suppl_internal_branches.pdf",
        fig_fld,
    )
    fig_internal_branches(
        b_dfs["internal"],
        "n. mutations",
        "n_muts",
        "suppl_internal_branches_muts.pdf",
        fig_fld,
    )

    rates = evaluate_rates(edf, b_dfs["internal"], aln_L, aln_L_core)
    fig_rates_overview(rates, fig_fld)
