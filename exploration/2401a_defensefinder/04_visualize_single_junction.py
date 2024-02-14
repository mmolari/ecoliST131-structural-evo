# %%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pathlib
import json
import pypangraph as pp
from Bio import Phylo

fig_fld = pathlib.Path("figs/f04")
fig_fld.mkdir(exist_ok=True, parents=True)

dset = "ST131_ABC"
fld = pathlib.Path(f"../../results/{dset}")

# load annotations
fld_data = fld / "annotations/junct_pos_asm20-100-5"
is_df = pd.read_csv(fld_data / "ISEScan_real.csv", index_col=0)
pp_df = pd.read_csv(fld_data / "genomad_real.csv", index_col=0)
in_df = pd.read_csv(fld_data / "integronfinder_real.csv", index_col=0)
df_df = pd.read_csv(fld_data / "defensefinder_real.csv", index_col=0)

# load info
df_info = pd.read_csv(
    fld / "annotations/defense_finder/systems_summary.tsv",
    sep="\t",
    index_col=0,
)


# load joint coordinates dictionary
pos_file = fld / "backbone_joints/asm20-100-5/joints_pos.json"
with open(pos_file) as f:
    jp = json.load(f)

# load junction info dataframes
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
df_j = pd.read_csv(df_file, index_col=0)
df_j["delta_len"] = df_j["max_length"] - df_j["min_length"]

# load tree
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
strains = [x.name for x in tree.get_terminals()]

# genome lengths
iso_L_file = fld / "pangraph/genome_lengths.csv"
iso_L = pd.read_csv(iso_L_file, index_col=0)["length"].to_dict()

df_j["n_defsys"] = (
    df_j.index.map(df_df["junction"].value_counts()).fillna(0).astype(int)
)
df_j["n_IS"] = df_j.index.map(is_df["junction"].value_counts()).fillna(0).astype(int)
df_j["n_prophages"] = (
    df_j.index.map(pp_df["junction"].value_counts()).fillna(0).astype(int)
)
df_j["n_integrons"] = (
    df_j.index.map(in_df["junction"].value_counts()).fillna(0).astype(int)
)

# %%
mask = df_j["n_iso"] == 222
mask = df_j["n_defsys"] > 0
# mask = df_j["n_categories"] == 2
df_j[mask].to_csv("data/junctions_defsys.csv")
df_j[mask]

# %%
junctions = []
for j in df_j[mask].index:
    junctions.append(j)

for junction in junctions:

    print(f"Junction: {junction}")
    print(df_j.loc[junction])
    # print("defense systems:")
    # mask = df_df["junction"] == junction
    # for idx, row in df_df[mask].iterrows():
    #     print(idx)
    #     print(row)

    j_pos = jp[junction]
    graph_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{junction}.json"
    pan = pp.Pangraph.load_json(graph_file)
    bdf = pan.to_blockstats_df()

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

        # plot ISs, prophages, integrons and defense systems
        el_displayed = []
        el_displayed.append((pp_df, "aqua"))
        el_displayed.append((in_df, "fuchsia"))
        el_displayed.append((is_df, "green"))
        el_displayed.append((df_df, "black"))
        for el_df, c in el_displayed:
            mask = el_df["junction"] == junction
            mask &= el_df["iso"] == iso

            for _, row in el_df[mask].iterrows():
                B, E = row["ib"], row["ie"]
                if E < B:
                    print(f"Iso: {iso} ie < ib {y=}")
                    continue
                if not strand:
                    E, B = -B, -E
                E, B = (E - delta) % L, (B - delta) % L
                if np.abs(E - L) < np.abs(E):
                    E -= L
                if np.abs(B - L) < np.abs(B):
                    B -= L

                yp = y
                axs[1].plot([B, E], [yp, yp], "|-", color=c, alpha=0.5)

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

    axs[1].set_title("annotations")
    axs[1].set_xlabel("position (bp)")
    axs[2].set_title("block structure")
    axs[2].set_xlabel("position (bp)")

    plt.tight_layout()
    plt.savefig(fig_fld / f"J_{junction}.png", dpi=200)
    plt.show()

    # save colors
    # block_color = {k: mpl.colors.to_hex(v) for k, v in block_color.items()}
    # cdf = pd.DataFrame.from_dict(block_color, orient="index", columns=["Color"])
    # cdf.to_csv(fig_fld / f"J_structure_{junction}_colors.csv")

# %%
# plot element legend
fig, ax = plt.subplots(1, 1, figsize=(5, 1))
names = ["prophage", "integron", "IS", "defense system"]
for i, (el_df, c) in enumerate(el_displayed):
    ax.plot([i, i + 0.4], [0, 0], "|-", color=c, lw=2, markersize=13)
    ax.text(i, 0.3, names[i], va="center")
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(-0.5, 0.5)
ax.axis("off")
plt.tight_layout()
plt.savefig(fig_fld / "element_legend.png", dpi=200)
plt.show()


# %%
