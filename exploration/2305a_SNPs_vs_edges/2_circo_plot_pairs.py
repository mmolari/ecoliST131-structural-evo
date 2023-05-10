# %%

import utils as ut

import numpy as np
from Bio import SeqIO, AlignIO, Phylo
from pycirclize import Circos


import pathlib

import matplotlib.pyplot as plt

fig_p = pathlib.Path("fig")
fig_p.mkdir(exist_ok=True)


def svfig(svname):
    for k in ["pdf", "png", "svg"]:
        plt.savefig(fig_p / f"{svname}.{k}", dpi=300, facecolor="w")


prefix = "../../results/ST131/pangraph"

aln_file = f"{prefix}/asm20-100-5-alignment/filtered_corealignment.fa"
tree_file = f"{prefix}/asm20-100-5-filtered-coretree.nwk"
pangraph_file = f"{prefix}/asm20-100-5-polished.json"

tree = Phylo.read(tree_file, format="newick")
tree.ladderize()
str_order = [l.name for l in tree.get_terminals()]

aln = AlignIO.read(aln_file, "fasta")
A = np.array(aln)
aln_strains = [s.name for s in aln]
# %%

# extract SNPs pairs
snps_pairs = ut.aln_to_branchings(aln_file)
snps_pairs = [b for b in snps_pairs if len(b) == 2]
# %%
block_pa = ut.extract_block_PA_df(pangraph_file, len_thr=100, remove_dupl=False)

# %%

bl_pairs = []
for bl, pa in block_pa.iteritems():
    vc = pa.value_counts()
    if len(vc) != 2:
        continue
    if vc.min() != 2:
        continue
    idxs = pa.index[pa == vc.index[1]]
    bl_pairs.append(list(idxs))
# %%

# extract edge dictionary
E_dict = ut.extract_edge_dictionary(pangraph_file, len_thr=100, remove_dupl=False)
strains_set = set(aln_strains)
N = len(strains_set)
# %%
edge_pairs = []
for k, s in E_dict.items():
    if len(s) == 2:
        edge_pairs.append(list(s))
    if len(s) == N - 2:
        edge_pairs.append(list(strains_set - s))


# %%


# %%

# Initialize circos sector with tree size
circos = Circos(sectors={"Tree": tree.count_terminals()})
sector = circos.sectors[0]

# Plot tree with inner style
track = sector.add_track((75, 100))
track.axis(ec="lightgrey")
track.tree(
    tree,
    outer=False,
    use_branch_length=True,
    leaf_label_size=0,
    line_kws=dict(color="black", lw=1),
    text_kws=dict(color="gray"),
)

circos.text(
    "\n".join(
        [
            f"n. pairwise",
            f"snps={len(snps_pairs)}",
            f"blocks={len(bl_pairs)}",
            f"edges={len(edge_pairs)}",
        ]
    ),
    size=10,
    deg=315,
    r=120,
)

for p1, p2 in snps_pairs:
    i1 = str_order.index(p1) + 0.5
    i2 = str_order.index(p2) + 0.5
    reg1 = ("Tree", i1, i1)
    reg2 = ("Tree", i2, i2)
    circos.link(reg1, reg2, color="red", lw=0.5, alpha=0.2, height_ratio=0.1)

for p1, p2 in bl_pairs:
    i1 = str_order.index(p1) + 0.5
    i2 = str_order.index(p2) + 0.5
    reg1 = ("Tree", i1, i1)
    reg2 = ("Tree", i2, i2)
    circos.link(reg1, reg2, color="blue", lw=0.5, alpha=0.05, height_ratio=0.2)

for p1, p2 in edge_pairs:
    i1 = str_order.index(p1) + 0.5
    i2 = str_order.index(p2) + 0.5
    reg1 = ("Tree", i1, i1)
    reg2 = ("Tree", i2, i2)
    circos.link(reg1, reg2, color="green", lw=0.5, alpha=0.03, height_ratio=0.23)


fig = circos.plotfig()
svfig("pairs_on_tree")

# %%
