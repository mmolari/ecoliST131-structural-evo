# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pypangraph as pp
from Bio import SeqIO
import pandas as pd
import pathlib

fig_fld = pathlib.Path("figs/f04b")
fig_fld.mkdir(parents=True, exist_ok=True)

fname = "../../results/ST131_ABC/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(fname)
bdf = pan.to_blockstats_df()
# ref = "NZ_CP096110.1"
ref = "NZ_CP124372.1"
len_thr = 500
# %%

fname = "../../exploration/2402b_figures_draft/figs/f04/ST131_ABC/block_stats.csv"
bc = pd.read_csv(fname, index_col=0)
fname = "../../exploration/2402b_figures_draft/figs/f04/ST131_ABC/mergers.csv"
mg = pd.read_csv(fname, index_col=0)
mg = mg["0"].to_dict()

# %%
fname = f"../../data/gbk/{ref}.gbk"
features = []
with open(fname, "r") as f:
    for record in SeqIO.parse(f, "genbank"):
        seq = record.seq
        for feat in record.features:
            if feat.type == "CDS":
                if "gene" in feat.qualifiers:
                    gene = feat.qualifiers["gene"][0]
                    if gene in ["dnaA", "tus"]:
                        features.append(feat)
                        print(feat)

# %%


def gc_skew(seq, bin_width):
    n = len(seq)
    bins = int(n / bin_width)
    skew = np.zeros(bins)
    for i in range(bins):
        start = i * bin_width
        end = (i + 1) * bin_width
        skew[i] = (seq[start:end].count("G") - seq[start:end].count("C")) / bin_width
    return skew


plt.plot(gc_skew(seq, 10000))
plt.axhline(0, color="gray", lw=0.5)
plt.show()

# %%


# extract paths and filter to core paths
def keep_f(bid):
    return bdf.loc[bid, "core"] and (bdf.loc[bid, "len"] > len_thr)


Bs = pan.paths[ref].block_ids
offset = pan.paths[ref].block_positions[0]
L_ref = sum([bdf["len"][b] for b in Bs])
# polar coordinates
fig, ax = plt.subplots(
    figsize=(6, 6),
    subplot_kw={"projection": "polar"},
)


factor = 2 * np.pi / L_ref
x = offset * factor
for bid in Bs:
    if not keep_f(bid):
        c = "gray"
        lw = 1
        # continue
    elif bid in mg:
        meta_block = mg[bid]
        c = bc.loc[meta_block, "color"]
        lw = 3
    else:
        c = bc.loc[bid, "color"]
        lw = 3
    l = bdf.loc[bid, "len"]
    # plt.plot([x, x + l], [0, 0], color=c, lw=lw)
    l *= factor
    plt.barh(y=10, width=l, height=lw, left=x, color=c, edgecolor="none")
    x += l

# add skew
bw = 20000
y_set = 7.5
y_skew = gc_skew(seq, bw)
y_skew = y_skew * 12 + y_set
x_skew = np.arange(len(y_skew), dtype=float) * bw
x_skew *= factor
# periodic boundary
y_skew = np.concatenate([y_skew, y_skew[:1]])
x_skew = np.concatenate([x_skew, x_skew[:1]])
ax.plot(x_skew, y_skew, color="gray", lw=1)
# densify
x_skew_d = np.linspace(0, 2 * np.pi, 1000)
y_skew_d = np.interp(x_skew_d, x_skew, y_skew, period=2 * np.pi)
ax.fill_between(
    x_skew_d, y_skew_d, y_set, color="C0", alpha=0.5, where=y_skew_d > y_set
)
ax.fill_between(
    x_skew_d, y_skew_d, y_set, color="C3", alpha=0.5, where=y_skew_d < y_set
)

# set y axis
ax.set_yticks([y_set])
ax.set_yticklabels([""])
ax.set_ylim([0, 10])

y_set = 6
for feat in features:
    loc = feat.location
    b = feat.location.start
    theta = b * factor
    gene = feat.qualifiers["gene"][0]
    ax.scatter([theta], [y_set], marker="o", color="gray")
    ax.text(theta, y_set - 1, gene, color="gray", ha="center", va="center")
xt = np.arange(0, L_ref, int(1e6))
ax.set_xticks(xt * factor)
ax.set_xticklabels([f"{x/1e6:.0f} Mb" for x in xt])
ax.set_title(f"Synteny blocks in {ref}")
plt.tight_layout()
plt.savefig(fig_fld / f"synteny_circle_{ref}.png", dpi=300)
plt.show()


# %%
