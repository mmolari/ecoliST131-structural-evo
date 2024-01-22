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

fig_fld = pathlib.Path("figs/f04")
fig_fld.mkdir(exist_ok=True, parents=True)


data_fld = pathlib.Path("data")
pp_to_j = pd.read_csv(data_fld / "phages_to_joints.csv", index_col=0)

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

df_j["n_prophages"] = pp_to_j["junction"].value_counts()
df_j["n_prophages"].fillna(0, inplace=True)
df_j["has_prophage"] = df_j["n_prophages"] > 0
df_j["prophage_freq"] = (
    pp_to_j.groupby("junction").apply(lambda x: len(np.unique(x["iso"])))
    / df_j["n_iso"]
)
df_j["prophage_freq"].fillna(0, inplace=True)

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
sns.scatterplot(data=df_j, x="nonempty_freq", y="prophage_freq")

mask = df_j["nonempty_freq"] < 0.1
# mask &= df_j["nonempty_freq"] < 0.1
mask &= df_j["prophage_freq"] < 0.1
# mask &= np.abs(df_j["prophage_freq"] - df_j["nonempty_freq"]) < 0.1

sns.scatterplot(data=df_j[mask], x="nonempty_freq", y="prophage_freq")

df_j[mask].sort_values("prophage_freq", ascending=False).head(20)

# %%
sns.scatterplot(data=df_j, x="delta_len", y="n_categories", alpha=0.1)
plt.xscale("log")
plt.yscale("log")

mask = df_j["prophage_freq"] > 0
# mask = df_j["nonempty_freq"] < 0.1
# mask &= df_j["nonempty_freq"] < 0.1
# mask &= np.abs(df_j["prophage_freq"] - df_j["nonempty_freq"]) < 0.1
sns.scatterplot(data=df_j[mask], x="delta_len", y="n_categories")
plt.show()


df_j[mask].sort_values("n_categories", ascending=False).head(20)
# %%

# junction = "GPXHVRRZLC_f__ZLAJFQLBFQ_r"
# junction = "CGQZOLYSTD_r__ZLLQQUXUIP_r"
# junction = "KAIPBNCIHR_r__WXCHSHHCDT_f"
# junction = "QJRASUKHLX_r__UHYGUNDBFL_f"
# junction = "CXKKDGMPSE_f__KBLPANZOCZ_f"
# junction = "SKPHAXSFLS_f__TORJAESOXF_f"
# junction = "IRXHOEIDDO_f__RVLRLPDSYQ_f"
# junction = "ATFHQYFPNW_f__GUDQDOMJFJ_f"
# junction = "XXVMWZCEKI_r__YUOECYBHUS_r"  # hotspot
# junction = "JVNRLCFAVD_f__PLTCZQCVRD_r"  # hotspot
# junction = "EJPOGALASQ_f__KUIFCLFQSI_r"  # hotspot
junction = "BWEZXGGFBK_r__MVMOFPVELT_r"  # hotspot

j_pos = jp[junction]

graph_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{junction}.json"
pan = pp.Pangraph.load_json(graph_file)
bdf = pan.to_blockstats_df()

mask = pp_to_j["junction"] == junction
pp_to_j[mask]
# %%

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
    if not iso in j_pos:
        continue
    Cb, Ab, Ae, Ce, strand = j_pos[iso]
    if Ce < Cb:
        print(f"Iso: {iso} Ce < Cb")
        continue

    if not strand:
        Cb, Ab, Ae, Ce = -Ce, -Ae, -Ab, -Cb
    delta = Cb
    Cb, Ab, Ae, Ce = Cb - delta, Ab - delta, Ae - delta, Ce - delta
    axs[1].plot([Cb, Ab], [y, y], lw=3, color="C0")
    axs[1].plot([Ab, Ae], [y, y], lw=3, color="lightgray")
    axs[1].plot([Ae, Ce], [y, y], lw=3, color="C1")

    # plot prophages
    mask = pp_to_j["junction"] == junction
    mask &= pp_to_j["iso"] == iso

    for _, row in pp_to_j[mask].iterrows():
        pb, pe = row["pb"], row["pe"]
        if pe < pb:
            print(f"Iso: {iso} pe < pb")
            continue
        if not strand:
            pe, pb = -pb, -pe
        pe, pb = pe - delta, pb - delta

        yp = y
        axs[1].plot([pb, pe], [yp, yp], "|-", color="red", alpha=0.5)

    # plot block structure
    path = pan.paths[iso]
    x = 0
    yb = y + 0.6
    for block_id in path.block_ids:
        l = bdf.loc[block_id, "len"]
        if block_id in block_color:
            c = block_color[block_id]
        else:
            c = cmap(len(block_color) % 20)
            block_color[block_id] = c

        axs[2].plot([x, x + l], [yb, yb], color=c, linewidth=3)
        x += l

axs[1].set_title("prophage annotations")
axs[1].set_xlabel("position (bp)")
axs[2].set_title("block structure")
axs[2].set_xlabel("position (bp)")

plt.tight_layout()
plt.savefig(fig_fld / f"J_structure_{junction}.png", dpi=200)
plt.show()
# %%
