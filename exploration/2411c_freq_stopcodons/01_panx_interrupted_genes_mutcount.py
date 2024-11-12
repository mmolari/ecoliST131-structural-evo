# %%
import pandas as pd
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
import utils as ut
import pypangraph as pp
from collections import Counter
from Bio import Seq, SeqIO, Phylo, AlignIO, Align
import subprocess

# fig_fld = pathlib.Path("figs/f01")
# fig_fld.mkdir(parents=True, exist_ok=True)

data_fld = pathlib.Path("res/f01")
data_fld.mkdir(parents=True, exist_ok=True)

tree_file = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
strains = [x.name for x in tree.get_terminals()]


gene_cluster_json = "../../data/panX/data/ST131_ABC/vis/geneCluster.json"
gc = ut.load_genecluster_json(gene_cluster_json)
gids = gc.index.to_numpy()

loc_fld_real = pathlib.Path("../../results/ST131_ABC/panx/gc_j_pos/asm20-100-5/real")


def parse_locations(gids, path):
    gc_info = []

    for gid in gids:
        res = {}
        row = gc.loc[gid]
        res["len"] = row["geneLen"]
        res["gid"] = gid
        res["n_iso"] = row["count"]
        res["n_loci"] = len(row["locus"].split())
        res["dupl"] = row["dupli"] == "yes"
        loc_file = path / f"genecl_{gid}.csv"
        if not loc_file.exists():
            # print(f"missing {loc_file}")
            raise ValueError

        try:
            df = pd.read_csv(loc_file)
            jc = df["junction"].value_counts()
            res["max_j"] = jc.max()
            res["j_nuniq"] = len(jc)

            res["acc"] = len(df)
            res["core"] = res["n_loci"] - len(df)
        except pd.errors.EmptyDataError:
            # file is empty
            # print(f"empty {loc_file}")
            res["core"] = res["n_loci"]

        res["ann"] = row["ann"]
        gc_info.append(res)
    gc_info = pd.DataFrame(gc_info)
    gc_info.set_index("gid", inplace=True)
    gc_info.fillna({"acc": 0, "max_j": 0, "j_nuniq": 0}, inplace=True)
    gc_info = gc_info.astype({"acc": int, "max_j": int, "j_nuniq": int})
    return gc_info


print("parsing real positions")
gci = parse_locations(gids, loc_fld_real)


# %%


# check number of perfectly placed almost-core genes

mask = ~gci["dupl"]
mask &= gci["j_nuniq"] == 0
mask &= gci["n_loci"] == 221
# mask &= gci["n_loci"] > (222 / 2)
sdf = gci[mask]
sdf.sort_values("n_loci", ascending=False)


# %%

# load pangraph
pan_file = "../../results/ST131_ABC/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pan_file)
locator = pp.Locator(pan)
# %%


def load_gid_pos(gid):
    fname = f"../../results/ST131_ABC/panx/gc_loci/genecl_{gid}.csv"
    df = pd.read_csv(fname, index_col=0)
    return df


def locate_genecluster(pdf, locator):
    loc_df = []
    for loc_id, row in pdf.iterrows():
        iso, beg, end = row["iso"], row["beg"], row["end"]
        blocks, pos, occ = locator.find_interval(iso, beg, end)
        loc_df.append(
            {
                "id": loc_id,
                "iso": iso,
                "blocks": blocks,
                "pos": pos,
                "occ": occ,
            }
        )
    loc_df = pd.DataFrame(loc_df).set_index("id")
    return loc_df


def unique_bl(loc_df):
    blocks = loc_df["blocks"]
    uniq = np.unique(blocks)
    return uniq


def aln_limits(A, ldf):
    iso_pos = ldf.set_index("iso")["pos"].apply(lambda x: x[0])
    begs = []
    ends = []
    for iso in iso_pos.index:
        a = A[iso]
        a = np.array(list(a))
        idxs = np.arange(len(a)) - np.cumsum(a == "-")
        pb, pe = iso_pos.loc[iso]
        # find corresponding idxs:
        b = np.where(idxs == pb)[0][0]
        e = np.where(idxs == pe)[0][-1]
        begs.append(b)
        ends.append(e)
    beg = min(begs)
    end = max(ends) - 1
    return beg, end


def stripped_aln(A, strains, b, e, ldf, complement):
    aln = []
    gene_strains = ldf["iso"].unique()
    for s in strains:
        a = list(A[s])[b:e]
        a = "".join(a)
        record = Seq.Seq(a)
        if complement:
            record = record.reverse_complement()
        d = "" if s in gene_strains else "missing"
        record = SeqIO.SeqRecord(record, id=s, description=d)
        aln.append(record)
    return aln


def decide_rev_comp(ldf):
    first_row = ldf.iloc[0]
    iso = first_row["iso"]
    locus_id = first_row.name
    block_strand = first_row["occ"][0][2]

    locus_df_f = f"../../data/loci_df/{iso}.csv"
    locus_strand = pd.read_csv(locus_df_f, index_col=0).loc[locus_id]["strand"]
    assert locus_strand in [1, -1]
    locus_strand = locus_strand == 1

    complement = block_strand != locus_strand
    return complement


def create_aln(pan, strains, block, ldf):
    beg = loc_df["pos"].apply(lambda x: x[0][0])
    end = loc_df["pos"].apply(lambda x: x[0][1])

    print(f"alignment for {block}")
    aln, occs = pan.blocks[block].alignment.generate_alignments()
    A = {o[0]: a for a, o in zip(aln, occs)}
    b, e = aln_limits(A, ldf)

    # reverse complement?
    complement = decide_rev_comp(ldf)
    print(f"complement: {complement}")

    A = stripped_aln(A, strains, b, e, ldf, complement)
    A = Align.MultipleSeqAlignment(A)
    return A


def stop_codon(A):
    sc_isos = []
    for a in A:
        iso = a.id

        # split in triplets removing gaps
        seq = str(a.seq)
        seq = seq.replace("-", "")
        triplets = [seq[i : i + 3] for i in range(0, len(seq), 3)]
        # check for stop codon
        stop = any([t in ["TAA", "TAG", "TGA"] for t in triplets[:-1]])
        if stop:
            sc_isos.append(iso)
    return sc_isos


def frameshift(A, threshold):
    aln_l = len(A[0])
    N = len(A)
    frames = np.zeros((N, aln_l), dtype=int)
    for n, a in enumerate(A):
        a = np.array(list(a.seq))
        frames[n] = np.cumsum(a == "-") % 3

    # consensus frame
    cons_frames = np.array(
        [Counter(frames[:, l]).most_common(1)[0][0] for l in range(aln_l)]
    )
    frames = frames != cons_frames
    # plt.matshow(frames)
    # plt.show()

    diffs = np.sum(frames, axis=1) / aln_l
    # plt.hist(diffs)
    # plt.show()

    fs_isos = []
    for n, d in enumerate(diffs):
        if d > threshold:
            fs_isos.append(A[n].id)
    return fs_isos


def analyze_aln(A, verbose=False):
    sc = stop_codon(A)
    fs = frameshift(A, threshold=0.1)
    missing = set([a.id for a in A if "missing" in a.description])

    sc_attributable = set(missing) == set(sc)
    fs_attributable = set(missing) == set(fs)

    if not fs_attributable and not sc_attributable:
        print(f"stop codon: {sc}")
        print(f"frameshift: {fs}")
        print(f"missing: {missing}")

    return {
        "stop_codon": sc_attributable,
        "frameshift": fs_attributable,
    }


res_df = []
for gid in sdf.index:
    print(f"\nprocessing {gid}")
    pdf = load_gid_pos(gid)

    loc_df = locate_genecluster(pdf, locator)
    uniq = unique_bl(loc_df)
    if len(uniq) > 1:
        print(f"## skipping, multiple paths {uniq}")
        continue
    block = uniq[0]
    if len(block) > 1:
        print(f"## skipping, multiple blocks {block}")
        continue
    block = block[0]

    # create and save alignment
    aln = create_aln(pan, strains, block, loc_df)
    f_save = data_fld / f"{gid}_aln.fa"
    SeqIO.write(aln, f_save, "fasta")

    # analyze alignment
    aln_res = analyze_aln(aln)

    res = {"gid": gid, "n_uniq": len(uniq), "uniq": uniq, "n_bl": len(uniq[0])}
    res |= aln_res
    res_df.append(res)

res_df = pd.DataFrame(res_df)
res_df[["n_uniq", "n_bl"]].value_counts()

# %%
res_df[["stop_codon", "frameshift"]].value_counts()

# %%
res_df.to_csv(data_fld / "stats.csv")
