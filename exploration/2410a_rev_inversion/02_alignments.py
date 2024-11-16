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
isolates = pan.strains()

iso_A = "NZ_CP107114.1"
iso_B = "NZ_CP107184.1"
iso_C = "NZ_CP107182.1"

# %%


def save_aln(M, isolates, fname):
    records = [
        SeqRecord.SeqRecord(Seq.Seq("".join(row)), id=iso, description="")
        for iso, row in zip(isolates, M)
    ]
    with open(fname, "w") as f:
        SeqIO.write(records, f, "fasta")


for mergers, prefix in [
    (["RYYAQMEJGY"], "merger_alignment"),
    (["TORJAESOXF", "PLTCZQCVRD"], "flanking_alignment"),
]:
    for merger in mergers:
        mg_df[merger] = merger

    # get the alignment for the merger
    path = pan.paths[iso_A]
    A = {iso: "" for iso in isolates}
    for bid, s in zip(path.block_ids, path.block_strands):
        if not bdf.loc[bid, "core"]:
            continue
        if bid not in mg_df.index:
            continue
        if mg_df.loc[bid] not in mergers:
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
    alignment = Align.MultipleSeqAlignment(records)

    # remove gaps
    M = np.array(alignment)
    mask = np.all(M != "-", axis=0)
    M = M[:, mask]

    lengths = {
        "L_full": len(mask),
        "L_ungapped": M.shape[1],
    }

    # export full alignment
    save_aln(M, isolates, res_fld / f"{prefix}_ungapped.fa")

    # remove all consensus sites
    mask = (M == M[0]).all(axis=0)
    M = M[:, ~mask]

    lengths["L_restr"] = M.shape[1]

    # export restricted alignment
    save_aln(M, isolates, res_fld / f"{prefix}_restricted.fa")

    ldf = pd.DataFrame(lengths, index=[0])
    ldf.to_csv(res_fld / f"{prefix}_lengths.csv", index=False)
# %%
