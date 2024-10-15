# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import pathlib
from collections import defaultdict
import pypangraph as pp
from Bio import SeqIO, AlignIO, SeqRecord, Align, Seq

# look into inversion that happened multiple times
# 115 kbps
# isolates:
# NZ_CP107114.1
# NZ_CP107184.1
# NZ_CP107182.1

fld = pathlib.Path("../../results/ST131_ABC")
res_fld = pathlib.Path("res")
res_fld.mkdir(exist_ok=True)

colors_file = fld / "pangraph/coresynt-asm20-100-5/blocks.csv"
mergers_file = fld / "pangraph/coresynt-asm20-100-5/mergers.csv"

cl_df = pd.read_csv(colors_file, index_col=0)
mg_df = pd.read_csv(mergers_file, index_col=0)["0"]

pangraph_file = fld / "pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()
# %%
iso_A = "NZ_CP107114.1"
iso_B = "NZ_CP107184.1"
iso_C = "NZ_CP107182.1"
merger = "RYYAQMEJGY"

mg_df["RYYAQMEJGY"] = "RYYAQMEJGY"
# %%
isolates = pan.strains()
# get the alignment for the merger
path = pan.paths[iso_A]
A = {iso: "" for iso in isolates}
for bid, s in zip(path.block_ids, path.block_strands):
    if not bdf.loc[bid, "core"]:
        continue
    if bid not in mg_df.index:
        continue
    if mg_df.loc[bid] != merger:
        continue
    block = pan.blocks[bid]
    S, O = block.alignment.generate_alignments()
    for seq, (iso, num, strand) in zip(S, O):
        if not s:
            seq = Seq.Seq(seq).reverse_complement()
            seq = str(seq)
        A[iso] += seq

# transform in biopython alignment
records = [SeqRecord.SeqRecord(Seq.Seq(A[iso]), id=iso) for iso in isolates]
# export
with open(res_fld / "merger_alignment.fa", "w") as f:
    SeqIO.write(records, f, "fasta")

# %%
# remove gaps
alignment = Align.MultipleSeqAlignment(records)
M = np.array(alignment)
mask = np.all(M != "-", axis=0)
M = M[:, mask]
M.shape

# export restricted alignment
records = [
    SeqRecord.SeqRecord(Seq.Seq("".join(row)), id=iso) for iso, row in zip(isolates, M)
]
with open(res_fld / "merger_alignment_no_gaps.fa", "w") as f:
    SeqIO.write(records, f, "fasta")
# %%
# %%
# %%
# %%
