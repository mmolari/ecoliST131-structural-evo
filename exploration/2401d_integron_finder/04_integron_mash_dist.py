# %%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pathlib
import seaborn as sns
import json
import subprocess
from Bio import SeqIO, Phylo
from collections import defaultdict

fig_fld = pathlib.Path("figs/f04")
fig_fld.mkdir(exist_ok=True, parents=True)


data_fld = pathlib.Path("data")
df_int = pd.read_csv(data_fld / "integrons_to_joints.csv", index_col=0)
df_int.set_index("integron", inplace=True)

dset = "ST131_ABC"

# load joint coordinates dictionary
fld = pathlib.Path(f"../../results/{dset}")
pos_file = fld / "backbone_joints/asm20-100-5/joints_pos.json"
with open(pos_file) as f:
    jp = json.load(f)

tree = Phylo.read(fld / "pangraph/asm20-100-5-filtered-coretree.nwk", "newick")
strains = [x.name for x in tree.get_terminals()]

# %%
records = []
for idx, row in df_int.iterrows():
    iso = row["iso"]
    ib, ie = row["ib"], row["ie"]
    fa_fname = fld / f"../../data/fa/{iso}.fa"
    with open(fa_fname) as f:
        seq = SeqIO.read(f, "fasta")
        seq = seq.seq[ib:ie]
    rec = SeqIO.SeqRecord(seq, id=idx, description="")
    records.append(rec)
mash_infile = data_fld / "integrons.fa"
with open(mash_infile, "w") as f:
    SeqIO.write(records, f, "fasta")

# %%
# mash distance
msh_file = data_fld / "integrons.msh"
cmd = f"mash triangle -i -k 21 -s 5000 {mash_infile} > {msh_file}"
subprocess.run(cmd, shell=True, check=True)


# %%
def parse_mash_file(fname):
    with open(fname) as f:
        lines = f.readlines()
    N = int(lines[0].strip())
    D = np.zeros((N, N))
    ids = []
    for n, l in enumerate(lines[1:]):
        L = l.strip().split("\t")
        ids.append(L[0])
        D[n, : len(L[1:])] = np.array(L[1:]).astype(float)
    D += D.T
    # create distance df
    D_df = pd.DataFrame(D, index=ids, columns=ids)
    return D_df


D_df = parse_mash_file(msh_file)
plt.matshow(D_df)
plt.colorbar()
plt.tight_layout()
plt.show()

# %%

order = df_int.sort_values(["type", "junction", "iso"], ascending=False).index


def color_generator_discr(palette, other_palette=None):
    cmap = mpl.colormaps[palette]
    for n in range(20):
        yield cmap(n)
    if other_palette is not None:
        cmap = mpl.colormaps[other_palette]
        for n in range(20):
            yield cmap(n)
    while True:
        yield "black"


tp_color = {"complete": "C0", "CALIN": "C1"}
J_cgen = color_generator_discr("tab10")
iso_cgen = color_generator_discr("tab20b", "tab20c")
J_color = defaultdict(lambda: next(J_cgen))
iso_color = defaultdict(lambda: next(iso_cgen))

fig, axs = plt.subplots(
    1, 2, figsize=(10, 7), sharey=True, gridspec_kw={"width_ratios": [1, 4]}
)
ax = axs[0]
for n, integron in enumerate(order):
    row = df_int.loc[integron]
    tp = row["type"]
    J = row["junction"]
    iso = row["iso"]
    ax.scatter(0, n + 0.5, color=tp_color[tp], s=100)
    ax.scatter(1, n + 0.5, color=J_color[J], s=100)
    ax.scatter(2, n + 0.5, color=iso_color[iso], s=100)

ax.set_xticks([0, 1, 2])
ax.set_xticklabels(["type", "junction", "iso"])
ax.set_xlim(-0.5, 2.5)

ax = axs[-1]
N = len(order)
new_ax = ax.inset_axes([0.8, 0.3, 0.02, 0.6])
sns.heatmap(
    D_df.loc[order, order],
    mask=np.triu(np.ones((N, N), dtype=bool), k=1),
    ax=ax,
    cmap="viridis_r",
    square=False,
    cbar_kws={"shrink": 0.5, "label": "Mash distance"},
    vmin=0,
    vmax=0.1,
    cbar_ax=new_ax,
)
ax.set_xticklabels([], rotation=90)


sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "integron_alignment.png", dpi=300)
plt.show()

# %%
tree = Phylo.read(fld / "pangraph/asm20-100-5-filtered-coretree.nwk", "newick")
# prune
for x in tree.get_terminals():
    if x.name not in iso_color:
        tree.prune(x)

fig, ax = plt.subplots(figsize=(3, 5))
label_func = lambda x: "" if x.name not in iso_color else x.name
color_func = lambda x: "black" if x not in iso_color else iso_color[x]
Phylo.draw(tree, do_show=False, label_func=label_func, label_colors=color_func, axes=ax)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
for spine in ax.spines.values():
    spine.set_visible(False)
plt.tight_layout()
plt.savefig(fig_fld / "integron_restricted_coretree.png", dpi=300)
plt.show()

# %%
