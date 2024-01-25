# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
import numpy as np
import json

import pypangraph as pp
import utils as ut

from collections import defaultdict
from Bio import Phylo, Seq, SeqIO, SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


dset = "ST131_ABC"
edge = "VCAVVOUNDI_f__XIWJABIXEM_f"
fig_fld = pathlib.Path(f"figs/n00b/{dset}/{edge}")
fig_fld.mkdir(exist_ok=True, parents=True)


def svfig(name):
    # plt.savefig(fig_fld / f"{name}.svg")
    plt.savefig(fig_fld / f"{name}.png", facecolor="white", dpi=150)


fld = pathlib.Path(f"../../results/{dset}")
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
df_edge_len = fld / "backbone_joints/asm20-100-5/edge_len.csv"
pos_file = fld / "backbone_joints/asm20-100-5/joints_pos.json"
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
gbk_fld = pathlib.Path(f"../../data/gbk/")

with open(pos_file) as f:
    jp = json.load(f)
df = pd.read_csv(df_file, index_col=0)
df_el = pd.read_csv(df_edge_len, index_col=0)
df["delta_len"] = df["max_length"] - df["min_length"]

tree = Phylo.read(tree_file, "newick")
isolates = [l.name for l in tree.get_terminals()]


# %%

gbk_files = {i: gbk_fld / f"{i}.gbk" for i in isolates}


def extract_sequence_and_genes(genbank_file, start_location, end_location):
    # Initialize variables to store results
    gene_annotations = []
    S, E = start_location, end_location

    loop_around = E < S

    # Iterate through each record in the GenBank file
    for record in SeqIO.parse(genbank_file, "genbank"):
        # Extract the sequence within the specified range (considering circular genomes)
        if loop_around:
            seq = record.seq[S - 1 :] + record.seq[:E]
        else:
            seq = record.seq[S - 1 : E]

        # Iterate through each feature in the record
        for feature in record.features:
            # Check if the feature is a gene and falls within or partially within the specified range
            if feature.type == "CDS":
                s = feature.location.nofuzzy_start
                e = feature.location.nofuzzy_end
                if s > e:
                    raise ValueError("s > e")
                if e - s > 1e6:
                    continue
                if (
                    ((s <= S) and (e > S))
                    or ((s < E) and (e >= E))
                    or (s >= S and e <= E)
                ):
                    # Extract the gene sequence and add it to the list
                    gene_annotations.append(
                        {
                            "locus_tag": feature.qualifiers["locus_tag"][0],
                            "product": feature.qualifiers["product"][0],
                            "strand": feature.location.strand,
                            "start": s,
                            "end": e,
                            "sub_start": s - S,
                            "sub_end": e - S,
                        }
                    )

        break
    return seq, pd.DataFrame(gene_annotations)


seqs = []
genes_df = []
for iso in isolates:
    pos = jp[edge][iso]
    print(iso)
    print(pos)
    sc, sa, ea, ec, strand = pos
    print(gbk_files[iso])
    flank = 500
    seq, genes = extract_sequence_and_genes(gbk_files[iso], sa - flank, ea + flank)
    if strand:
        seq = seq.reverse_complement()

    seqs.append(SeqRecord.SeqRecord(seq, id=iso, description=""))
    genes["iso"] = iso
    genes_df.append(genes)

SeqIO.write(seqs, fig_fld / f"{edge}.fa", "fasta")
genes_df = pd.concat(genes_df)
genes_df.to_csv(fig_fld / f"{edge}_genes.csv")


# %%

# run mafft
# mafft --thread 8 --auto --keeplength --addfragments {edge}.fa {edge}.fa > {edge}_mafft.fa
import subprocess

in_file = fig_fld / f"{edge}.fa"
out_file = fig_fld / f"{edge}.aln.fa"
subprocess.run(
    f"mafft --thread 4 --auto {in_file} > {out_file}",
    shell=True,
)

# %%
# open with aliview
subprocess.run(f"aliview {out_file}", shell=True)

# %%
genes_df["product"].value_counts()
# %%
import subprocess

pg_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{edge}.json"
fg_file = fig_fld / "subgraph.png"
subprocess.run(
    f"""
    python3 plot_junction_categories.py \
            --pangraph {pg_file} \
            --tree {tree_file} \
            --fig {fg_file} \
            --len_filt 0
    """,
    shell=True,
)

# %%
