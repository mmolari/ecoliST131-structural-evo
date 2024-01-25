# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
import numpy as np
import json

import pypangraph as pp
import utils as ut

from Bio import Phylo, SeqIO


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


dset = "ST131_ABC"
fig_fld = pathlib.Path(f"figs/n1/{dset}")
fig_fld.mkdir(exist_ok=True, parents=True)


def svfig(name):
    # plt.savefig(fig_fld / f"{name}.svg")
    plt.savefig(fig_fld / f"{name}.png", facecolor="white")


fld = pathlib.Path(f"../../results/{dset}")
df_file = fld / "backbone_joints/asm20-100-5/junctions_stats.csv"
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
joints_pos_file = fld / "backbone_joints/asm20-100-5/joints_pos.json"

# %%
df = pd.read_csv(df_file, index_col=0)

mask = df["n_iso"] > 50
mask &= df["n_categories"] == 2
df = df[mask]
sdf = df[df["singleton"]].copy()

# %%

tree = Phylo.read(tree_file, "newick")
tree.root_at_midpoint()
tree.ladderize()
terminal_len = {b.name: b.branch_length for b in tree.get_terminals()}
isolates = [t.name for t in tree.get_terminals()]

with open(f"data/{dset}/singleton_gainloss_info.json", "r") as f:
    gainloss_info = json.load(f)

edf = pd.read_csv(f"data/{dset}/singleton_event_df.csv")

with open(joints_pos_file, "r") as f:
    joints_pos = json.load(f)

joints_pos = {k: joints_pos[k] for k in gainloss_info.keys()}

# %%
ann_df = []
for j, gli in gainloss_info.items():
    iso = gli["iso"]
    event_type = gli["type"]
    iso = gli["iso"]
    block_names = gli["event_blocks"]
    print(f"---\n{iso=} | {event_type=} | {block_names=}")

    # only gain/loss can be determined
    if event_type == "other":
        continue
    elif event_type == "loss":
        # for losses, need to change the reference iso
        print(f"loss: {iso=}")
        I = set([i for i in isolates if i in joints_pos[j]])
        iso = np.random.choice(list(I - set(iso)))
        print(f"random iso: {iso=}")

    # unpack joint position
    core_start, start_acc, end_acc, core_end, strand = joints_pos[j][iso]
    print(f"{j=} {iso=}")
    print(f"{core_start=} {start_acc=} {end_acc=} {core_end=} {strand=}")

    # load block name and info

    # load pangraph
    pan_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{j}.json"
    pan = pp.Pangraph.load_json(pan_file)

    # gbk file
    gbk_file = f"../../data/gbk/{iso}.gbk"
    gbk_L = len(list(SeqIO.parse(gbk_file, "genbank"))[0].seq)

    for b in block_names:
        path = pan.paths[iso]
        bids = path.block_ids
        where_b = np.where(bids == b)[0].flatten()
        assert path.offset is None
        for i in where_b:
            beg, end = path.block_positions[i : i + 2]
            print(f"{beg=} {end=}")
            if strand:
                beg += core_start
                end += core_start
            else:
                beg, end = (core_end - end, core_end - beg)

            beg, end = beg % gbk_L, end % gbk_L

            print(f"extracting {beg=} {end=}")
            if beg > end:
                print("WARNING: beg > end")
                _, genes1 = extract_sequence_and_genes(gbk_file, beg, gbk_L)
                _, genes2 = extract_sequence_and_genes(gbk_file, 0, end)
                genes = pd.concat([genes1, genes2])
            else:
                _, genes = extract_sequence_and_genes(gbk_file, beg, end)
            if len(genes) == 0:
                continue
            genes["block"] = b
            genes["block_len"] = len(pan.blocks[b])
            genes["event_type"] = event_type
            genes["iso"] = iso
            genes["junction"] = j
            L = genes["sub_end"] - genes["sub_start"]
            out = -np.minimum(genes["sub_start"], 0)
            out += np.maximum(genes["sub_end"] - genes["block_len"], 0)
            # print(out)
            # print(L)
            genes = genes[(out / L) < 0.5]
            ann_df.append(genes)

ann_df = pd.concat(ann_df)
ann_df = ann_df.reset_index(drop=True)
ann_df.to_csv(f"data/{dset}/singleton_junctions_annotations.csv", index=False)


# %%

rnd_df = []
for j, gli in gainloss_info.items():
    iso = gli["iso"]
    event_type = gli["type"]
    iso = gli["iso"]
    block_names = gli["event_blocks"]
    print(f"---\n{iso=} | {event_type=} | {block_names=}")

    # only gain/loss can be determined
    if event_type == "other":
        continue
    elif event_type == "loss":
        # for losses, need to change the reference iso
        print(f"loss: {iso=}")
        I = set([i for i in isolates if i in joints_pos[j]])
        iso = np.random.choice(list(I - set(iso)))
        print(f"random iso: {iso=}")

    # unpack joint position
    core_start, start_acc, end_acc, core_end, strand = joints_pos[j][iso]
    print(f"{j=} {iso=}")
    print(f"{core_start=} {start_acc=} {end_acc=} {core_end=} {strand=}")

    # load block name and info

    # load pangraph
    pan_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{j}.json"
    pan = pp.Pangraph.load_json(pan_file)

    # gbk file
    gbk_file = f"../../data/gbk/{iso}.gbk"
    gbk_L = len(list(SeqIO.parse(gbk_file, "genbank"))[0].seq)

    for b in block_names:
        path = pan.paths[iso]
        bids = path.block_ids
        where_b = np.where(bids == b)[0].flatten()
        assert path.offset is None
        for i in where_b:
            beg, end = path.block_positions[i : i + 2]
            print(f"{beg=} {end=}")
            if strand:
                beg += core_start
                end += core_start
            else:
                beg, end = (core_end - end, core_end - beg)

            rand_delta = np.random.randint(gbk_L)
            beg += rand_delta
            end += rand_delta
            beg, end = beg % gbk_L, end % gbk_L

            print(f"extracting {beg=} {end=}")
            if beg > end:
                print("WARNING: beg > end")
                _, genes1 = extract_sequence_and_genes(gbk_file, beg, gbk_L)
                _, genes2 = extract_sequence_and_genes(gbk_file, 0, end)
                genes = pd.concat([genes1, genes2])
            else:
                _, genes = extract_sequence_and_genes(gbk_file, beg, end)
            if len(genes) == 0:
                continue
            genes["block"] = b
            genes["block_len"] = len(pan.blocks[b])
            genes["event_type"] = event_type
            genes["iso"] = iso
            genes["junction"] = j
            L = genes["sub_end"] - genes["sub_start"]
            out = -np.minimum(genes["sub_start"], 0)
            out += np.maximum(genes["sub_end"] - genes["block_len"], 0)
            # print(out)
            # print(L)
            genes = genes[(out / L) < 0.5]
            rnd_df.append(genes)

rnd_df = pd.concat(rnd_df)
rnd_df = rnd_df.reset_index(drop=True)
rnd_df.to_csv(f"data/{dset}/singleton_junctions_annotations_random.csv", index=False)

# %%

ann_df = pd.read_csv(f"data/{dset}/singleton_junctions_annotations.csv")
rnd_df = pd.read_csv(f"data/{dset}/singleton_junctions_annotations_random.csv")
# %%

gain_df = ann_df[ann_df["event_type"] == "gain"].copy()
loss_df = ann_df[ann_df["event_type"] == "loss"]

# check whether gains have transposases or phage proteins
pr = gain_df["product"].str.lower()
gain_df["h_tr"] = (
    pr.str.contains("transposase")
    | pr.str.contains("transposon")
    | pr.str.contains("insertion sequence")
)
gain_df["h_ph"] = pr.str.contains("phage")
gain_df["h_mges"] = gain_df["h_tr"] | gain_df["h_ph"]

# for every junction, check if any transposon, IS or phage is present
ggj = gain_df.groupby("junction")
j_df = {
    "transp.": ggj["h_tr"].sum() > 0,
    "prophage": ggj["h_ph"].sum() > 0,
    "MGE": ggj["h_mges"].sum() > 0,
}
j_df = pd.DataFrame(j_df)

rnd_gain_df = rnd_df[rnd_df["event_type"] == "gain"].copy()

# check whether gains have transposases or phage proteins
pr = rnd_gain_df["product"].str.lower()
rnd_gain_df["h_tr"] = (
    pr.str.contains("transposase")
    | pr.str.contains("transposon")
    | pr.str.contains("insertion sequence")
)
rnd_gain_df["h_ph"] = pr.str.contains("phage")
rnd_gain_df["h_mges"] = gain_df["h_tr"] | gain_df["h_ph"]

# for every junction, check if any transposon, IS or phage is present
ggj = rnd_gain_df.groupby("junction")
rj_df = {
    "transp.": ggj["h_tr"].sum() > 0,
    "prophage": ggj["h_ph"].sum() > 0,
    "MGE": ggj["h_mges"].sum() > 0,
}
rj_df = pd.DataFrame(rj_df)

# %%
N = len(list(filter(lambda x: x["type"] == "gain", gainloss_info.values())))

X = {"data": j_df.sum() / N, "reshuffled data": rj_df.mean() / N}
X = pd.DataFrame(X).T
X = X.unstack().reset_index()
X.columns = ["element type", "source", "fraction of joints with element"]

fig, ax = plt.subplots(1, 1, figsize=(6, 4))
sns.barplot(data=X, x="element type", y="fraction of joints with element", hue="source", ax=ax)
plt.yscale("log")
sns.despine()
plt.tight_layout()
svfig("mge_annotations_log")
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(6, 4))
sns.barplot(data=X, x="element type", y="fraction of joints with element", hue="source", ax=ax)
sns.despine()
plt.tight_layout()
svfig("mge_annotations")
plt.show()


# %%
j_df["len"] = gain_df.groupby("junction")["block_len"].mean()
# %%
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
sns.histplot(
    data=j_df, x="len", hue="prophage", multiple="stack", log_scale=True, ax=ax
)
plt.xlabel("junction length")
plt.title("gain length distribution")
sns.despine()
plt.tight_layout()
svfig("prophage_vs_len")
plt.show()

# %%
