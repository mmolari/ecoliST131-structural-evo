# %%
# - [ ] check from blast sequence right features
# - [ ] check sum of length of features vs length of window. Outliers?
# - [ ] check IS or transposases per window. (integrase/transposase/IS/recombinase...)
# - [ ] check betalactam
# - [ ] check phages

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

df_ann = pd.read_csv("data/annotations.csv", index_col=0)
df_ev = pd.read_csv("data/twocat_events.csv", index_col=0)

# %%
df_ann.columns
# %%
L1 = df_ann.groupby("joint")["feature_len"].sum()
L2 = df_ev["block len"]
# %%
x = pd.DataFrame([L1, L2]).T
# %%
sns.scatterplot(data=x, x="feature_len", y="block len")
plt.xscale("log")
plt.yscale("log")
plt.show()
# %%

# check for beta-lactam
df_ann["CTX-M-15"] = df_ann["product"].str.lower().str.contains("ctx-m-15")
df_ann["betalactam"] = df_ann["product"].str.lower().str.contains("beta-lactam")
# check for MGE
terms = ["transposase", "integrase", "recombinase", "conjug"]
df_ann["mge"] = df_ann["product"].str.lower().str.contains("|".join(terms), regex=True)
df_ann["phage"] = df_ann["product"].str.lower().str.contains("phage")
# merge
for k in ["betalactam", "CTX-M-15", "mge", "phage"]:
    df_ev[k] = df_ann.groupby("joint")[k].any()

df_ev["mge_fract"] = df_ann.groupby("joint")["mge"].mean()

# %%
df_ev[["betalactam", "CTX-M-15", "mge", "phage"]].value_counts()
# %%
gbk_fld = "../../data/gbk"
strains = pd.read_csv("data/mugration/PA_states.csv", index_col=0).index.to_list()

from Bio import SeqIO
from collections import defaultdict

mge_fract = defaultdict(list)

print(f"n. strains = {len(strains)}")
for iso in strains:
    gbk_file = f"{gbk_fld}/{iso}.gbk"
    print(f"reading {gbk_file}")
    for rec in SeqIO.parse(gbk_file, "genbank"):
        print(f"record {rec.id} {rec.description}")
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers["product"]
            q = q[0].lower()
            has_mge = any([q.find(t) != -1 for t in terms])
            mge_fract[iso].append(has_mge)
        break
mge_fract = {k: np.mean(v) for k, v in mge_fract.items()}
mge_df = pd.DataFrame.from_dict({"mge_fraction": mge_fract})
mge_df.to_csv("data/mge_fract.csv")

# %%
mge_df = pd.read_csv("data/mge_fract.csv")

sns.histplot(data=mge_df, x="mge_fraction", bins=20)
plt.axvline(mge_df["mge_fraction"].mean(), color="red", label="mean")
plt.xlabel("fraction of MGE annotations")
plt.legend()
sns.despine()
plt.tight_layout()
plt.savefig("figs/mge_fract_gbk.png", dpi=300)
plt.show()

# mge-fraction: between 1.5 and 2.5 percent

# %%
fig, axs = plt.subplots(2, 2, figsize=(10, 8))
ax = axs[0, 0]
sns.histplot(data=df_ev, x="type", hue="event category", ax=ax, multiple="stack")

# how many are associated with single blocks? (gain/loss) Only 6 "other"

ax = axs[0, 1]
sns.histplot(
    data=df_ev[df_ev["event category"] != "other"],
    x="block len",
    ax=ax,
    bins=15,
    log_scale=True,
    hue="event category",
    common_bins=True,
    stat="count",
    multiple="stack",
    palette=["C1", "C2"],
)

ax = axs[1, 0]
sns.histplot(data=df_ev, x="type", hue=df_ev["isolate"].isna(), multiple="stack", ax=ax)
ax.legend(["multiple events", "single event"])

ax = axs[1, 1]
sns.histplot(
    data=df_ev, x="event category", hue="mge", ax=ax, multiple="stack", palette="Set2"
)
ax.legend(["MGE", "no MGE"])

sns.despine()
plt.tight_layout()
plt.savefig("figs/annotation_stats.pdf")
plt.savefig("figs/annotation_stats.png")
plt.show()

# %%
from Bio import Phylo

tree_file = "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
tree.ladderize()
branch_len = {b.name: b.branch_length for b in tree.get_terminals()}
branch_len |= {b.name: b.branch_length for b in tree.get_nonterminals()}
dfbl = pd.DataFrame.from_dict(branch_len, orient="index", columns=["branch len"])

dfbl["n. events"] = df_ev.groupby("isolate").count().iloc[:, 0]
dfbl["n. events"] = dfbl["n. events"].fillna(0)
# %%

sns.scatterplot(data=dfbl, x="branch len", y="n. events", alpha=0.5)
plt.xscale("log")
plt.xlim(1e-7, 1e-4)
sns.despine()
plt.tight_layout()
plt.savefig("figs/annotation_stats.pdf")
plt.savefig("figs/annotation_stats.png")
plt.show()

# %%
