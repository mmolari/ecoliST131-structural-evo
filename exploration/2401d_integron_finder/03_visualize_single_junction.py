# %%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pathlib
import seaborn as sns
import json
import pypangraph as pp
from Bio import SeqIO, Phylo

fig_fld = pathlib.Path("figs/f03")
fig_fld.mkdir(exist_ok=True, parents=True)


data_fld = pathlib.Path("data")
pp_to_j = pd.read_csv(data_fld / "integrons_to_joints.csv", index_col=0)

dset = "ST131_ABC"

# load joint coordinates dictionary
fld = pathlib.Path(f"../../results/{dset}")
pos_file = fld / "backbone_joints/asm20-100-5/joints_pos.json"
with open(pos_file) as f:
    jp = json.load(f)

# load junction info dataframes
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
df_j = pd.read_csv(df_file, index_col=0)
df_j["delta_len"] = df_j["max_length"] - df_j["min_length"]

df_j["n_integrons"] = pp_to_j["junction"].value_counts()
df_j["n_integrons"].fillna(0, inplace=True)
df_j["has_integron"] = df_j["n_integrons"] > 0
df_j["integron_freq"] = (
    pp_to_j.groupby("junction").apply(lambda x: len(np.unique(x["iso"])))
    / df_j["n_iso"]
)
df_j["integron_freq"].fillna(0, inplace=True)

# load tree
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
strains = [x.name for x in tree.get_terminals()]


# isolates genome length
iso_L = {}
for iso in strains:
    fa_fname = fld / f"../../data/fa/{iso}.fa"
    with open(fa_fname) as f:
        iso_L[iso] = len(SeqIO.read(f, "fasta"))

# %%
sns.scatterplot(data=df_j, x="delta_len", y="n_categories", alpha=0.1)
plt.xscale("log")
plt.yscale("log")

mask = df_j["integron_freq"] > 0
# mask = df_j["nonempty_freq"] < 0.1
# mask &= df_j["nonempty_freq"] < 0.1
# mask &= np.abs(df_j["integron_freq"] - df_j["nonempty_freq"]) < 0.1
sns.scatterplot(data=df_j[mask], x="delta_len", y="n_categories")
plt.show()


df_j[mask].sort_values("n_categories", ascending=False).head(20)
# %%

# junction = "KGJXRPELGX_r__KIVDDRSTJR_r"  #   24
# junction = "UXPSNWWHML_f__YHLTNASHXN_f"  #    2
# junction = "DSLXSJZSTL_r__KGJXRPELGX_f"  #    2
# junction = "DSLXSJZSTL_r__QQSILILDBT_r"  #    2
# junction = "RYYAQMEJGY_r__ZTHKZYHPIX_f"  #    1
# junction = "KIVDDRSTJR_f__WDXIDBBCPP_r"  #    1
# junction = "HKKMVAAVLK_r__KJFGXXNSDH_r"  #    1
# junction = "KGJXRPELGX_r__NHYWGMABYN_f"  #    1
# junction = "OWCPAELREV_f__URGBYGHEXB_f"  #    1
# junction = "KGJXRPELGX_r__TCWDRAKLPS_r"  #    1
# junction = "EQHTXGCHGZ_r__KGJXRPELGX_f"  #    1

Junctions = [
    "KGJXRPELGX_r__KIVDDRSTJR_r",
    "UXPSNWWHML_f__YHLTNASHXN_f",
    "DSLXSJZSTL_r__KGJXRPELGX_f",
    "DSLXSJZSTL_r__QQSILILDBT_r",
    "RYYAQMEJGY_r__ZTHKZYHPIX_f",
]

for junction in Junctions:
    j_pos = jp[junction]

    graph_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{junction}.json"
    pan = pp.Pangraph.load_json(graph_file)
    bdf = pan.to_blockstats_df()

    mask = pp_to_j["junction"] == junction
    pp_to_j[mask]

    fig, axs = plt.subplots(
        1, 3, figsize=(15, 15), gridspec_kw={"width_ratios": [1, 3, 3]}, sharey=True
    )
    Phylo.draw(tree, axes=axs[0], do_show=False, label_func=lambda x: "")
    axs[0].set_title("Core tree")

    # plot
    block_color = {}
    cmap = mpl.colormaps.get_cmap("tab20")
    y = 0
    for iso in strains:
        y += 1
        L = iso_L[iso]
        if not iso in j_pos:
            continue
        Cb, Ab, Ae, Ce, strand = j_pos[iso]
        if Ce < Cb:
            print(f"Iso: {iso} Ce < Cb, {y=}")

        if not strand:
            Cb, Ab, Ae, Ce = -Ce, -Ae, -Ab, -Cb
        delta = Cb
        Cb, Ab, Ae, Ce = (
            (Cb - delta) % L,
            (Ab - delta) % L,
            (Ae - delta) % L,
            (Ce - delta) % L,
        )
        axs[1].plot([Cb, Ab], [y, y], lw=3, color="C0")
        axs[1].plot([Ab, Ae], [y, y], lw=3, color="lightgray")
        axs[1].plot([Ae, Ce], [y, y], lw=3, color="C1")

        # plot integrons
        mask = pp_to_j["junction"] == junction
        mask &= pp_to_j["iso"] == iso

        for _, row in pp_to_j[mask].iterrows():
            pb, pe = row["ib"], row["ie"]
            if pe < pb:
                print(f"Iso: {iso} pe < pb")
                continue
            if not strand:
                pe, pb = -pb, -pe
            pe, pb = (pe - delta) % L, (pb - delta) % L
            if np.abs(pe - L) < np.abs(pe):
                pe -= L
            if np.abs(pb - L) < np.abs(pb):
                pb -= L

            yp = y
            axs[1].plot([pb, pe], [yp, yp], "|-", color="red", alpha=0.5)

        # plot block structure
        path = pan.paths[iso]
        x = 0
        yb = y
        for block_id in path.block_ids:
            l = bdf.loc[block_id, "len"]
            if block_id in block_color:
                c = block_color[block_id]
            else:
                c = cmap(len(block_color) % 20)
                block_color[block_id] = c

            axs[2].plot([x, x + l], [yb, yb], color=c, linewidth=3)
            x += l

    axs[1].set_title("integron annotations")
    axs[1].set_xlabel("position (bp)")
    axs[2].set_title("block structure")
    axs[2].set_xlabel("position (bp)")

    plt.tight_layout()
    plt.savefig(fig_fld / f"J_structure_{junction}.png", dpi=200)
    plt.show()

    # save colors
    block_color = {k: mpl.colors.to_hex(v) for k, v in block_color.items()}
    cdf = pd.DataFrame.from_dict(block_color, orient="index", columns=["Color"])
    cdf.to_csv(fig_fld / f"J_structure_{junction}_colors.csv")
# %%
