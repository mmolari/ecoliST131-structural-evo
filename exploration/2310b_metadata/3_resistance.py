# %%
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from Bio import Phylo

tree_file = "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk"
res_file = "../../results/ST131/resistance/ncbi_summary.txt"

tree = Phylo.read(tree_file, "newick")
tree.ladderize()

# %%


def parse_element(x, thr):
    """
    Assign `.` = 0
    and `100.00;100;00;100.00` = 3
    """
    if x == ".":
        return 0
    x = str(x)
    els = x.split(";")
    ct = 0
    for el in els:
        if float(el) > thr:
            ct += 1
    return ct


def load_res_df(fname, thr):
    df = pd.read_csv(fname, sep="\t", index_col=0)
    # transform filename to accession number
    df.index = df.index.map(lambda x: pathlib.Path(x).stem)
    df.drop(columns="NUM_FOUND", inplace=True)
    # parse elements
    df = df.applymap(lambda x: parse_element(x, thr))
    return df


df = load_res_df(res_file, 95.0)
df
# %%
Y, X = df.shape

iso_y = {l.name: n + 1 for n, l in enumerate(tree.get_terminals())}

fig, axs = plt.subplots(
    1,
    2,
    figsize=(0.2 * X + 3, 1 + 0.08 * Y),
    sharey=True,
    gridspec_kw={"width_ratios": [1, X * 0.10]},
)

# plot tree
ax = axs[0]
Phylo.draw(tree, lambda x: "", do_show=False, axes=ax)

# get colormap
cmap = mpl.cm.get_cmap("rainbow")
norm = mpl.colors.Normalize(vmin=1, vmax=df.max().max())
mapp = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)


# plot count matrix
ax = axs[1]
xticks = {}
x = 0
for col in df.columns:
    xticks[x] = col
    for iso in tree.get_terminals():
        y = iso_y[iso.name]
        ct = df.loc[iso.name, col]
        if ct > 0:
            c = cmap(norm(ct))
            ax.scatter(x, y, color=c, marker="s", edgecolor="black")
    ax.axvline(x, c="lightgray", zorder=-1, alpha=0.3)
    x += 1
ax.set_xticks(list(xticks.keys()))
ax.set_xticklabels(list(xticks.values()), rotation=90)
ax.set_title("n. resistance genes")
plt.colorbar(mapp, ax=ax, shrink=0.1)


for ax in axs:
    for n in range(Y):
        ax.axhline(n + 1, c="lightgray", zorder=-1, alpha=0.3)
    ax.grid(True, alpha=0.3, axis="y")
plt.tight_layout()
plt.savefig("figs2/resistance.png", dpi=300, facecolor="white")
plt.show()

# %%
