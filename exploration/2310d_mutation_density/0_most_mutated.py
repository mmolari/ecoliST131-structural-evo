# %%

from Bio import AlignIO, SeqIO, Phylo
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt

dset = "ST131_full"
subfld = f"../../results/{dset}/pangraph"
tree_file = f"{subfld}/asm20-100-5-filtered-coretree.nwk"
alignment_file = f"{subfld}/asm20-100-5-alignment/corealignment.fa"
filtered_alignment_file = f"{subfld}/asm20-100-5-alignment/filtered_corealignment.fa"
aln_info = f"{subfld}/asm20-100-5-alignment/corealignment_info.json"
filt_aln_info = f"{subfld}/asm20-100-5-alignment/filtered_corealignment_info_size.json"


def consensus(A):
    """given an alignment matrix returns the consensus"""

    def site_consensus(l):
        lett, ct = np.unique(l, return_counts=True, axis=0)
        return lett[np.argmax(ct)]

    return np.apply_along_axis(site_consensus, 0, A)


def n_consensus_muts(aln):
    A = np.array(aln)

    cons = consensus(A)

    N_muts = (A != cons).sum(axis=1)
    N_muts = {a.id: n for a, n in zip(aln, N_muts)}
    return N_muts


aln = AlignIO.read(alignment_file, "fasta")
N_muts = n_consensus_muts(aln)
with open(aln_info, "r") as f:
    data = json.load(f)
    aln_L = data["n. consensus"] + data["n. snps"]

aln_f = AlignIO.read(filtered_alignment_file, "fasta")
N_muts_f = n_consensus_muts(aln_f)
with open(filt_aln_info, "r") as f:
    data = json.load(f)
    aln_f_L = data["polished aln size"]
# %%
tree = Phylo.read(tree_file, "newick")
tree.root_at_midpoint()
tree.ladderize()
# %%

leaves = [l.name for l in tree.get_terminals()]
L = len(leaves)

fig, axs = plt.subplots(1, 3, figsize=(8, L * 0.1), sharey=True)

ax = axs[0]
Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")

ax = axs[1]
ax.barh(np.arange(L) + 1, [N_muts[l] / aln_L for l in leaves])
ax.set_xlabel("distance to consensus")
ax.set_title("raw alignment")

ax = axs[2]
ax.barh(np.arange(L) + 1, [N_muts_f[l] / aln_f_L for l in leaves])
ax.set_xlabel("distance to consensus")
ax.set_title("filtered alignment")

plt.tight_layout()
plt.savefig("figs/consensus_mutations.png", dpi=300, facecolor="w")
plt.show()
# %%
