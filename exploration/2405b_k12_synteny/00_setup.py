# %%
from Bio import SeqIO, SeqRecord, Seq
import subprocess
import pypangraph as pp

ref = "NZ_CP096110.1"
thr_len = 500
# %%

src = f"../../data/fa/{ref}.fa"
dst = f"data/genomes/{ref}.fa"
cmd = ["cp", src, dst]
subprocess.run(cmd)

# %%
fname = "../../results/ST131_ABC/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(fname)
bdf = pan.to_blockstats_df()

# find core blocks
mask = bdf["core"] & (bdf["len"] >= 500)
CB = bdf[mask].index
# %%

records = []
for cb in CB:
    block = pan.blocks[cb]
    seq = Seq.Seq(block.sequence)
    rec = SeqRecord.SeqRecord(seq, id=cb, description="")
    records.append(rec)

with open("data/coreblocks.fa", "w") as f:
    SeqIO.write(records, f, "fasta")


# %%

# ============ blast search =============

# make databases
mk_cmd = lambda i, o: [
    ". $CONDA_PREFIX/../../etc/profile.d/conda.sh &&" "conda activate blast &&",
    "makeblastdb",
    f"-in {i}",
    "-dbtype nucl",
    f"-out {o}",
]

cmd1 = mk_cmd("data/genomes/k12.fa", "data/blastdb/k12/k12")
subprocess.run(" ".join(cmd1), shell=True)

cmd2 = mk_cmd(f"data/genomes/{ref}.fa", f"data/blastdb/{ref}/{ref}")
subprocess.run(" ".join(cmd2), shell=True)
# %%

# blast map
mk_cmd = lambda db, q, o: [
    ". $CONDA_PREFIX/../../etc/profile.d/conda.sh &&",
    "conda activate blast &&",
    "blastn",
    # "-subject_besthit",
    f"-db {db}",
    f"-query {q}",
    '-outfmt "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen slen"',
    f"-out {o}",
]
cmd1 = mk_cmd(
    "data/blastdb/k12/k12",
    "data/coreblocks.fa",
    "data/blast_res/k12.tsv",
)
subprocess.run(" ".join(cmd1), shell=True)

cmd2 = mk_cmd(
    f"data/blastdb/{ref}/{ref}",
    "data/coreblocks.fa",
    f"data/blast_res/{ref}.tsv",
)
subprocess.run(" ".join(cmd2), shell=True)

cmd3 = mk_cmd(
    f"data/blastdb/{ref}/{ref}",
    "data/genomes/k12.fa",
    "data/blast_res/k12_on_ref.tsv",
)
subprocess.run(" ".join(cmd3), shell=True)

# %%
