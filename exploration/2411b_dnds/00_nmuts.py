# %%
import json
import pandas as pd
import numpy as np
from Bio import AlignIO, Seq
import gzip
from functools import cache
from collections import Counter
import pathlib

res_fld = pathlib.Path("res")
res_fld.mkdir(exist_ok=True)


# load gene clusters info
gene_cluster_json = "../../data/panX/data/ST131_ABC/vis/geneCluster.json"

with open(gene_cluster_json, "r") as f:
    gdf = pd.DataFrame(json.load(f))
gdf["divers"] = gdf["divers"].astype(float)
gdf["count"] = gdf["count"].astype(int)
gdf["geneLen"] = gdf["geneLen"].astype(int)
gdf["event"] = gdf["event"].astype(int)
gdf.set_index("geneId", inplace=True)

# %%
N = 222
core_mask = gdf["count"] == 222
core_mask &= gdf["dupli"] == "no"
gdf["core"] = core_mask
gdf[core_mask]

# 2653 core gene clusters

# %%
# parse alignments


def dnds_alignment_count(aln):
    # turn into a matrix
    M = np.array(aln)
    assert M.shape[1] % 3 == 0, "Alignment length not a multiple of 3"

    N_aa = M.shape[1] // 3

    # process each codon
    bck_ct, mts_ct, tot_ct = Counter(), Counter(), Counter()
    for i in range(N_aa):
        codon = M[:, i * 3 : (i + 1) * 3]
        res = process_codon_aln(codon)
        if res is None:
            continue
        bck, tot, muts = res
        bck_ct += bck
        mts_ct += muts
        tot_ct += tot

    return bck_ct, tot_ct, mts_ct


@cache
def codon_to_aa(cod):
    cod = "".join(cod)
    aa = Seq.Seq(cod).translate()
    return str(aa)


def codon_consensus(C):
    # check if there are any gaps
    if "-" in C:
        print("gap")
        return None

    # take the most common nucleotide per position
    consensus = np.array([""] * 3, dtype=str)
    for i in range(3):
        col = C[:, i]
        vals, counts = np.unique(col, return_counts=True)
        idx = np.argmax(counts)
        consensus[i] = vals[idx]
    return tuple(consensus)


def n_polymorphic_columns(C):
    npol = 0
    for i in range(3):
        col = C[:, i]
        vals, counts = np.unique(col, return_counts=True)
        if len(vals) > 1:
            npol += 1
    return npol


def analyze_mutations(C, nt_cons):
    keep = True

    aa_cons = codon_to_aa(nt_cons)

    tot_ct = Counter()
    tot_ct["ncol"] = 3
    tot_ct["npol"] = n_polymorphic_columns(C)

    mut_ct = Counter()

    if tot_ct["npol"] == 0:
        return tot_ct, mut_ct, keep

    vals, cts = np.unique(C, return_counts=True, axis=0)
    for v in vals:
        v = tuple(v)
        if v != nt_cons:
            # is it synonymous?
            aa = codon_to_aa(v)
            syn = aa != aa_cons
            syn = "S" if syn else "N"
            tot_ct[syn] += 1
            # which nucleotide is changing
            nmuts = 0
            for i in range(3):
                if v[i] != nt_cons[i]:
                    nmuts += 1
                    mut_ct[(nt_cons[i], v[i], syn)] += 1
                    # break
            if nmuts > 1:
                print(f"More than one mutation, {v} vs {nt_cons}, discarding codon")
                keep = False

    return tot_ct, mut_ct, keep


@cache
def background_expectation(cons_nt):
    """count the effect of all possible mutations on a codon.
    Returns a counter whose keys are (nt, new_nt, syn)."""
    mut_counts = Counter()

    cons_aa = codon_to_aa(cons_nt)
    nts = set(list("ACGT"))
    for i in range(3):
        nt = cons_nt[i]
        for new_nt in nts - {nt}:
            # create mutated new codon
            new_cod = list(cons_nt)
            new_cod[i] = new_nt
            new_cod = tuple(new_cod)
            new_aa = codon_to_aa(new_cod)
            syn = cons_aa == new_aa
            syn = "S" if syn else "N"
            mut_counts[(nt, new_nt, syn)] += 1
    return mut_counts


def process_codon_aln(C):
    # check if there are any gaps
    if "-" in C:
        return None

    # get consensus and aa
    cons_nt = codon_consensus(C)

    bck = background_expectation(cons_nt)
    tot, muts, keep = analyze_mutations(C, cons_nt)

    if not keep:
        return None

    return bck, tot, muts


mt_counts = {}
bck_counts = {}
tot_counts = {}
for idx, row in gdf[core_mask].iterrows():
    gene_id = idx
    msa_idx = row["msa"]
    aln_file = f"../../data/panX/data/ST131_ABC/vis/geneCluster/{msa_idx}_na_aln.fa.gz"
    with gzip.open(aln_file, "rt") as f:
        aln = AlignIO.read(f, "fasta")
    bck_ct, tot_ct, mts_ct = dnds_alignment_count(aln)
    mt_counts[gene_id] = mts_ct
    bck_counts[gene_id] = bck_ct
    tot_counts[gene_id] = tot_ct

mt_counts = pd.DataFrame(mt_counts).fillna(0).astype(int)
bck_counts = pd.DataFrame(bck_counts).fillna(0).astype(int)
tot_counts = pd.DataFrame(tot_counts).fillna(0).astype(int)

# %%

mt_counts.to_csv(res_fld / "mt_counts.csv")
bck_counts.to_csv(res_fld / "bck_counts.csv")
tot_counts.to_csv(res_fld / "tot_counts.csv")

# %%
