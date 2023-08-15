# %%

from Bio import SeqIO
import pandas as pd
import json

import pypangraph as pp
import pathlib

# %%

genomes_fld = pathlib.Path("../../data/fa")
edge_pos_file = pathlib.Path(
    "../../results/ST131/backbone_joints/asm20-100-5/joints_pos.json"
)

pan = pp.Pangraph.load_json("data/pangraph.json")

# %%
with open(edge_pos_file, "r") as f:
    edge_pos = json.load(f)

# %%
for e, P in edge_pos.items():
    records = []
    Lmax = 0
    for iso, pos in P.items():
        beg, end, strand = pos
        with open(genomes_fld / f"{iso}.fa", "r") as f:
            seq = SeqIO.read(f, "fasta")
        if beg < end:
            sseq = seq[beg - 1 : end]
        else:
            sseq = seq[beg - 1 :] + seq[:end]
        if not strand:
            sseq = sseq.reverse_complement()
        sseq.id = f"{iso}"
        sseq.description = ""
        records.append(sseq)
        Lmax = max(Lmax, len(sseq))
    print(f"{e}: {Lmax/1e6}")

    # with open("test.fa", "w") as f:
    #     SeqIO.write(records, f, "fasta")
    # break
# %%

b = pan.blocks["CDCBSJDHBF"]
b_seq = b.sequence

with open("ref.fa", "w") as f:
    f.write(f">{b.id}\n{b_seq}\n")
# %%
