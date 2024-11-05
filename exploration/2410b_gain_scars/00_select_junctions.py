# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import pypangraph as pp
from Bio import SeqIO, AlignIO, Phylo
import subprocess
from collections import Counter

fld = pathlib.Path("../../results/ST131_ABC")
res_fld = pathlib.Path("res/j_cutouts")
res_fld.mkdir(parents=True, exist_ok=True)


terminal_df_fname = fld / "rates/asm20-100-5/terminal_coldspot.csv"
tdf = pd.read_csv(terminal_df_fname, index_col=0)

tree_fname = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_fname, "newick")
iso_order = [cl.name for cl in tree.get_terminals()]

# %%
mask = tdf["n_blocks"] == 3
mask &= tdf["isescan"] > 0
mask &= tdf["genomad"] == 0
mask &= tdf["event_type"] == "gain"

tdf[mask]

# %%


def get_junction_cutpoints(pan, isolates):
    paths = pan.to_paths_dict()

    first_block = [p[0] for p in paths.values()]
    last_block = [p[-1] for p in paths.values()]
    assert len(set(first_block)) == 1
    assert len(set(last_block)) == 1
    first_block = first_block[0]
    last_block = last_block[0]

    j_cutpoints = {}
    for iso in isolates:
        path = pan.paths[iso]
        pos = path.block_positions
        j_cutpoints[iso] = (pos[1], pos[-2])
    return j_cutpoints


def junctions_cutout(jcp, fname, delta):
    with open(fname, "r") as f:
        seqs = list(SeqIO.parse(f, "fasta"))
    cut_seqs = {}
    for seq in seqs:
        iso = seq.id
        start, end = jcp[iso]
        cut_seqs[iso] = seq.seq[start - delta : end + delta]
    cut_seqrecords = [
        SeqIO.SeqRecord(seq, id=iso, description="") for iso, seq in cut_seqs.items()
    ]
    return cut_seqrecords


for edge, row in tdf[mask].iterrows():
    print(edge)
    gain_iso = row["event_iso"]

    graph_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{edge}.json"
    seq_file = fld / f"backbone_joints/asm20-100-5/joints_seqs/{edge}.fa"
    pan = pp.Pangraph.load_json(graph_file)

    jcp = get_junction_cutpoints(pan, iso_order)
    csr = junctions_cutout(jcp, seq_file, 100)

    with open(res_fld / f"{edge}.fa", "w") as f:
        SeqIO.write(csr, f, "fasta")

    # execute alignment with mafft
    out_file = res_fld / f"{edge}.aln"
    cmd = f"conda run -n pangraph mafft --auto {res_fld / f'{edge}.fa'} > {out_file}"
    subprocess.run(cmd, shell=True)


# %%

# compare output files
res_df = []
for edge, row in tdf[mask].iterrows():
    gain_iso = row["event_iso"]

    # compare output files
    aln_file = res_fld / f"{edge}.aln"
    with open(aln_file, "r") as f:
        aln = list(AlignIO.read(f, "fasta"))

    # index of gain_iso
    gain_idx = [seq.id for seq in aln].index(gain_iso)
    A = np.array(aln)
    # remove gain_iso
    A = np.delete(A, gain_idx, axis=0)

    # look at how many different categories of sequences we have
    seq_dict = {seq.id: seq.seq for seq in aln}
    del seq_dict[gain_iso]

    # count the number of different sequences
    seq_counts = Counter(seq_dict.values())

    # count gaps per position
    gaps = np.sum(A == "-", axis=0)
    n_noncons_gaps = np.sum((gaps > 0) & (gaps < len(A)))

    res_df.append(
        {
            "edge": edge,
            "n_diff_seqs": len(seq_counts),
            "n_noncons_gaps": n_noncons_gaps,
        }
    )

res_df = pd.DataFrame(res_df)
res_df

# %%
mask = res_df["n_noncons_gaps"] > 0
res_df[mask]

# %%
# check the alignment
edges = res_df[mask]["edge"].values
for edge in edges:
    aln_file = res_fld / f"{edge}.aln"
    # open with aliview
    cmd = f"aliview {aln_file} && sleep 20"
    subprocess.run(cmd, shell=True)

# %%
