# %%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pathlib
import json
import pypangraph as pp
from Bio import Phylo

fig_fld = pathlib.Path("figs/f03")
fig_fld.mkdir(exist_ok=True, parents=True)

dset = "ST131_ABC"
fld = pathlib.Path(f"../../results/{dset}")

# load annotations
fld_data = fld / "annotations/junct_pos_asm20-100-5"
is_df = pd.read_csv(fld_data / "ISEScan_real.csv", index_col=0)
pp_df = pd.read_csv(fld_data / "genomad_real.csv", index_col=0)

# load info
is_info = pd.read_csv(fld / "annotations/isescan/is_summary.tsv", sep="\t")
is_info["id"] = (
    is_info["seqID"] + "|" + is_info["family"] + "|" + is_info.index.astype(str)
)
is_info = is_info.set_index("id", verify_integrity=True)


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

# %%
mask = df_j["n_categories"] == 2
mask &= df_j["n_iso"] == 222
df_j[mask].head(20)

# %%
junctions = []
junctions.append("KAIPBNCIHR_r__WXCHSHHCDT_f")  # prophage + hotspot
junctions.append("EJPOGALASQ_f__KUIFCLFQSI_r")  # hotspot

# # most IS junctions:
junctions.append("XXVMWZCEKI_r__YUOECYBHUS_r")  # N = 3125
junctions.append("RYYAQMEJGY_r__ZTHKZYHPIX_f")  # N = 2333
junctions.append("IFRPFFEGON_r__TCWDRAKLPS_r")  # N = 1476
junctions.append("JVNRLCFAVD_f__PLTCZQCVRD_r")  # N =  826
junctions.append("XFUKGTZLAV_f__YBWQKVQGZE_f")  # N =  684

# two categories
junctions.append("JKRDVEYDGL_f__OTPJRRWJRK_f")
junctions.append("JMOMDSHCBS_r__SPDPCYMYDN_r")  # independent 2 IS
junctions.append(
    "JHJVGMLICW_f__QHIDYNXVPV_f"
)  # nice one, exactly IS prediction! (no gene interruption)
junctions.append(
    "DSTCDJCESN_f__SBTAELODZT_f"
)  # nice one, exactly IS prediction! (no gene interruption)
junctions.append(
    "IXLMXEMXWI_r__XRXZJDDTTM_r"
)  # nice one, in subclade. Gene disruption: flagellar assembly peptidoglycan hydrolase FlgJ
junctions.append(
    "IZBYPELCGF_r__RVYEEDGLQI_r"
)  # nice one, longer IS, (GENE INTERRUPTED - inverse autotransporter invasin YchO, virulence factor?)


# IS and prophages
junctions.append("QNVVEFTBVQ_r__TLVFRBMGBC_r")  # 7 pp and 8 IS
junctions.append("GPXHVRRZLC_f__ZLAJFQLBFQ_r")  # 36 pp and 36 IS
junctions.append("CTWLDYAKJU_r__FQGXEWGAMU_r")  # 5 pp and 5 IS
junctions.append("DDVKABGVWS_r__IXLMXEMXWI_r")  # 1 pp and 1 IS (complete IS)
junctions.append("KSXZRIJFEH_r__SIHIJVBWQQ_r")  # 1 pp and 1 IS (complete IS)
junctions.append("SBTAELODZT_f__YYFLULYKQO_f")  # 1 pp and 1 IS (complete IS)


for junction in junctions:
    j_pos = jp[junction]
    graph_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{junction}.json"
    pan = pp.Pangraph.load_json(graph_file)
    bdf = pan.to_blockstats_df()

    mask = is_df["junction"] == junction
    print(is_df[mask])

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

        # plot ISs and prophages
        el_displayed = []
        el_displayed.append((pp_df, "black"))
        el_displayed.append((is_df, "red"))
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

    axs[1].set_title("IS annotations")
    axs[1].set_xlabel("position (bp)")
    axs[2].set_title("block structure")
    axs[2].set_xlabel("position (bp)")

    plt.tight_layout()
    plt.savefig(fig_fld / f"J_structure_{junction}_combined.png", dpi=200)
    plt.show()

    # save colors
    # block_color = {k: mpl.colors.to_hex(v) for k, v in block_color.items()}
    # cdf = pd.DataFrame.from_dict(block_color, orient="index", columns=["Color"])
    # cdf.to_csv(fig_fld / f"J_structure_{junction}_colors.csv")

# %%
