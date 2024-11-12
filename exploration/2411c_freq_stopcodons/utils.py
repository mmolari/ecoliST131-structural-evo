import pandas as pd
import numpy as np
from Bio import Seq, SeqIO, AlignIO
from functools import cache
from collections import Counter
import json


def mut_type(aa_old: str, aa_new: str) -> str:
    syn = aa_old == aa_new
    nonsense = aa_new == "*"
    return {
        (True, False): "S",
        (False, False): "M",
        (False, True): "*",
        (True, True): "?",
    }.get((syn, nonsense))


@cache
def codon_to_aa(cod: tuple) -> str:
    cod = "".join(cod)
    aa = Seq.Seq(cod).translate()
    return str(aa)


def n_polymorphic_columns(C: np.ndarray) -> int:
    """Returns the number of polymorphic columns in a codon sub-matrix."""
    npol = 0
    for i in range(3):
        col = C[:, i]
        vals = np.unique(col)
        if len(vals) > 1:
            npol += 1
    return npol


def codon_consensus(C: np.ndarray) -> tuple:
    """
    Given a codon sub-matrix, return the consensus codon as a tuple.
    Also check for gaps.
    """
    # check if there are any gaps
    has_gap = "-" in C

    # take the most common nucleotide per position
    consensus = np.array([""] * 3, dtype=str)
    for i in range(3):
        col = C[:, i]
        vals, counts = np.unique(col, return_counts=True)
        idx = np.argmax(counts)
        consensus[i] = vals[idx]
    return tuple(consensus), has_gap


@cache
def background_expectation(cons_nt: tuple) -> Counter:
    """Counts the effect of all possible mutations in a codon.
    Returns a Counter with the counts of each type of mutation.
    """
    possible_muts = Counter()

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
            mt = mut_type(cons_aa, new_aa)
            possible_muts[(nt, new_nt, mt)] += 1
    return possible_muts


def count_mutations(C, nt_cons):
    aa_cons = codon_to_aa(nt_cons)

    res = {
        "ncol": 3,
        "npol": n_polymorphic_columns(C),
        "muts": Counter(),
        "keep": True,
    }

    if res["npol"] == 0:
        return res

    vals, cts = np.unique(C, return_counts=True, axis=0)
    for v in vals:
        v = tuple(v)
        if v != nt_cons:
            # is it synonymous?
            aa = codon_to_aa(v)
            mt = mut_type(aa_cons, aa)

            # which nucleotide is changing
            nmuts = 0
            for i in range(3):
                if v[i] != nt_cons[i]:
                    nmuts += 1
                    res["muts"][(nt_cons[i], v[i], mt)] += 1
                    # break
            if nmuts > 1:
                print(f"More than one mutation, {v} vs {nt_cons}, {aa_cons} -> {aa}")
                # res["keep"] = False

    return res


def process_one_codon_aln(C):
    # get consensus and aa
    cons_nt, has_gaps = codon_consensus(C)

    # check if there are any gaps
    if has_gaps:
        return {"keep": False, "reason": "has_gaps"}

    bck = background_expectation(cons_nt)

    ct_muts = count_mutations(C, cons_nt)

    if not ct_muts["keep"]:
        return {"keep": False, "reason": "too_many_muts"}

    res = {
        "mut_possible": bck,
        "mut_count": ct_muts["muts"],
        "aln_stats": Counter({"npol": ct_muts["npol"], "ncol": ct_muts["ncol"]}),
        "keep": True,
    }
    return res


def process_alignment(aln: AlignIO.MultipleSeqAlignment) -> dict:
    """Takes as input a biopython alignment object and returns the counts of
    background mutations, total mutations, and total codons.
    """
    # turn into a matrix
    M = np.array(aln)
    assert M.shape[1] % 3 == 0, "Alignment length not a multiple of 3"

    N_aa = M.shape[1] // 3

    # process each codon
    res = {
        "mut_possible": Counter(),
        "mut_count": Counter(),
        "aln_stats": Counter(),
    }
    for i in range(N_aa):
        codon = M[:, i * 3 : (i + 1) * 3]
        codon_res = process_one_codon_aln(codon)
        if not codon_res["keep"]:
            match codon_res["reason"]:
                case "has_gaps":
                    res["aln_stats"]["n_gap_codons"] += 1
                case "too_many_muts":
                    res["aln_stats"]["n_too_many_muts"] += 1
                case _:
                    ValueError(f"Unknown reason: {codon_res['reason']}")
        else:
            res["mut_possible"] += codon_res["mut_possible"]
            res["mut_count"] += codon_res["mut_count"]
            res["aln_stats"] += codon_res["aln_stats"]

    return res


def gbk_to_loci_df(gbk_file):
    gb = SeqIO.read(gbk_file, "genbank")

    iso = gb.annotations["source"]

    # create a dataframe of loci
    loci = []
    for f in gb.features:
        if f.type == "CDS":
            loc = f.location
            # check if compount positions
            if len(loc.parts) > 1:
                # print(f.qualifiers["locus_tag"][0], "compound")
                # print(loc.parts)
                start = loc.parts[0].start
                end = loc.parts[-1].end
            else:
                start = loc.start
                end = loc.end
            loci.append(
                {
                    "locus": f.qualifiers["locus_tag"][0],
                    "start": start,
                    "end": end,
                    "strand": loc.strand,
                    "length": len(f),
                }
            )
    loci = pd.DataFrame(loci)
    return loci


def load_genecluster_json(fname):
    with open(fname, "r") as f:
        jc = pd.DataFrame(json.load(f))
    jc["geneId"] = jc["geneId"].astype(int)
    jc["divers"] = jc["divers"].astype(float)
    jc["count"] = jc["count"].astype(int)
    jc["geneLen"] = jc["geneLen"].astype(int)
    jc["event"] = jc["event"].astype(int)
    jc.set_index("geneId", inplace=True)
    return jc
