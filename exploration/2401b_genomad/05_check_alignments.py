# %%
import pandas as pd
import pathlib
from Bio import SeqIO, Phylo
import subprocess


dset = "ST131_ABC"


fig_fld = pathlib.Path("figs/f05")
fig_fld.mkdir(exist_ok=True, parents=True)


data_fld = pathlib.Path("data")
pp_to_j = pd.read_csv(data_fld / "phages_to_joints.csv", index_col=0)

dset = "ST131_ABC"

# load joint coordinates dictionary
fld = pathlib.Path(f"../../results/{dset}")


# load tree
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
strains = [x.name for x in tree.get_terminals()]


# isolates genome length
iso_L = {}
for iso in strains:
    fa_fname = fld / f"../../data/fa/{iso}.fa"
    with open(fa_fname) as f:
        iso_L[iso] = len(SeqIO.read(f, "fasta"))

data_fld = pathlib.Path("data/phage_aln")
data_fld.mkdir(exist_ok=True, parents=True)

# %%
# junction = "ANHYQUDQAA_r__IMAMHFLCWS_f"
junction = "ATFHQYFPNW_f__GUDQDOMJFJ_f"

mask = pp_to_j["junction"] == junction

records = {}
for _, row in pp_to_j[mask].iterrows():
    iso = row["iso"]
    pb, pe = row["pb"], row["pe"]
    fa_fname = fld / f"../../data/fa/{iso}.fa"
    with open(fa_fname) as f:
        seq = SeqIO.read(f, "fasta")
    phage_seq = seq.seq[pb:pe]
    if not row["strand"]:
        phage_seq = phage_seq.reverse_complement()
    print(f"{iso} {row['strand']} {pb} {pe} {len(phage_seq)}")
    if iso in records:
        raise ValueError(f"Duplicate iso {iso}")
    records[iso] = SeqIO.SeqRecord(
        phage_seq, id=iso, description=f"{row['strand']} {pb} {pe}"
    )

records = [records[iso] for iso in strains if iso in records]
# save as fasta file
with open(data_fld / f"{junction}.fa", "w") as f:
    SeqIO.write(records, f, "fasta")

# %%

out_file = data_fld / f"{junction}.aln.fa"
if not out_file.exists():
    # align using mafft
    command = f"mafft --auto --thread 8 {data_fld / junction}.fa > {out_file}"
    subprocess.run(command, shell=True, check=True)

# visualize with aliview
command = f"aliview {out_file}"
subprocess.run(command, shell=True, check=True)

# %%
