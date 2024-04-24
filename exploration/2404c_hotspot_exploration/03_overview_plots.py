# %%
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import pathlib

fig_fld = pathlib.Path("figs/f03")
fig_fld.mkdir(parents=True, exist_ok=True)

fname = "../../results/ST131_ABC/hotspots/asm20-100-5/hotspots.csv"
hdf = pd.read_csv(fname, index_col=0)
hdf
# %%
hsp_fld = pathlib.Path("../../results/ST131_ABC/hotspots/asm20-100-5/hotspot_stats")

# bins = np.logspace(-5, -3.7, 10)
bins = np.logspace(-6, -3.7, 10)
bins = [0] + list(bins)

vals = [
    "merge_n_shared",
    "merge_acc_len",
    "merge_core_len",
    "merge_n_edges",
    "n. breakpoints",
    "n. shared blocks",
    "n. private blocks",
]
for val in vals:
    fig, axs = plt.subplots(
        2,
        2,
        figsize=(8, 8),
        sharex="col",
        gridspec_kw={"width_ratios": [1, 0.05], "height_ratios": [0.3, 0.7]},
    )

    norm = mpl.colors.LogNorm(vmin=5e4, vmax=1e6)
    cmap = plt.cm.rainbow

    for hs, row in hdf.iterrows():
        # if row["pangenome_len"] < 1e5:
        #     continue
        # if row["n_categories"] < 20:
        #     continue

        fname = hsp_fld / f"{hs}.csv"
        df = pd.read_csv(fname, index_col=[0, 1])

        df["cut"] = pd.cut(df["core_div_filtered"], bins)
        gb = df.groupby("cut", observed=False)

        ax = axs[1, 0]
        x = gb["core_div_filtered"].mean()
        y = gb[val].mean()
        err = gb[val].std()
        col_val = row["pangenome_len"]
        col = cmap(norm(col_val))
        ax.errorbar(
            x,
            y,
            # yerr=err,
            marker="o",
            color=col,
            alpha=0.5,
            label=hs,
            linestyle="--",
            capsize=5,
        )

    # set axis off
    ax = axs[0, 1]
    ax.axis("off")

    # barplot count
    ax = axs[0, 0]
    ax.set_yscale("log")
    x = gb["core_div_filtered"].mean()
    n = gb["core_div_filtered"].count()
    ax.plot(x, n, "ko--")
    ax.set_ylabel("n. pairs")

    for b in bins:
        for ax in axs[:, 0]:
            ax.axvline(b, color="k", alpha=0.1)

    mapp = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapp.set_array([])
    cbar = plt.colorbar(mapp, cax=axs[1, 1])
    cbar.set_label("Pangenome length")
    ax = axs[1, 0]
    # ax.set_yscale("symlog")
    ax.set_xlabel("core genome divergence")
    ax.set_ylabel(val)
    # ax.set_xscale("symlog", linthresh=1e-5)
    ax.set_xscale("symlog", linthresh=1e-6)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    val_lab = val.replace(". ", "_").replace(" ", "_")
    plt.savefig(fig_fld / f"v2_{val}.png")
    plt.show()
# %%
