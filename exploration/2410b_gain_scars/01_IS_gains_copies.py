# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import pypangraph as pp
from Bio import SeqIO, AlignIO, Phylo
from collections import Counter

fld = pathlib.Path("../../results/ST131_ABC")
res_fld = pathlib.Path("res/01_IS_copies")
res_fld.mkdir(parents=True, exist_ok=True)


terminal_df_fname = fld / "rates/asm20-100-5/terminal_coldspot.csv"
tdf = pd.read_csv(terminal_df_fname, index_col=0)

tree_fname = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_fname, "newick")
iso_order = [cl.name for cl in tree.get_terminals()]

is_location_file = fld / "annotations/junct_pos_asm20-100-5/ISEScan_real.csv"
is_loc = pd.read_csv(is_location_file, index_col=0)

is_all_file = fld / "annotations/isescan/is_summary.tsv"
is_all = pd.read_csv(is_all_file, sep="\t")
is_all["n_other_copies"] = is_all["ncopy4is"] - 1

lengths_filename = fld / "pangraph/genome_lengths.csv"
Ls = pd.read_csv(lengths_filename, index_col=0)

# %%
mask = tdf["n_blocks"] == 3
mask &= tdf["isescan"] > 0
mask &= tdf["genomad"] == 0
mask &= tdf["event_type"] == "gain"

gains_df = tdf[mask]
# %%

mask = is_loc["junction"].isin(gains_df.index)
is_sub = is_loc[mask]

# %%

#  find the IS element

length_series = is_sub["iso"].map(Ls["length"])
ds = np.abs(is_sub["ib"] - is_sub["jab"])
de = np.abs(is_sub["ie"] - is_sub["jae"])
ds = np.minimum(ds, length_series - ds)
de = np.minimum(de, length_series - de)
delta = ds + de

plt.figure(figsize=(4, 3))
plt.hist(
    delta,
    bins=list(range(0, 300, 1)) + [delta.max() + 1],
    cumulative=True,
    histtype="step",
)

plt.xlim(0, 100)
plt.ylabel("n. ISs (cumulative)")
plt.xlabel("IS/junction distance (bp)")
sns.despine()
plt.tight_layout()
plt.savefig(res_fld / "IS_junction_distance.png", dpi=300)
plt.show()

# %%
event_iso = is_sub["junction"].map(gains_df["event_iso"])
is_sub[event_iso == is_sub["iso"]]
# %%

mask = delta < 40
mask &= event_iso == is_sub["iso"]
selected_ISs = is_sub[mask].index
is_sub[mask]
# %%
# export dataframe of selected ISs.

is_all["idx"] = (
    is_all["seqID"] + "|" + is_all["family"] + "|" + is_all.index.astype(str)
)
IS_gains = is_all.set_index("idx").loc[selected_ISs]

IS_gains
# %%
IS_gains["type"].value_counts()

IS_gains["n_other_copies"].value_counts()


# %%
def copies_histplot(df, ax):
    sns.histplot(
        data=df,
        x="n_other_copies",
        discrete=True,
        color="gray",
        element="step",
        ax=ax,
    )
    sns.histplot(
        data=df,
        x="n_other_copies",
        discrete=True,
        hue="family",
        hue_order=["IS1", "IS3", "IS66", "IS21"],
        multiple="stack",
        ax=ax,
    )
    ax.set_xticks(range(0, df["n_other_copies"].max() + 1, 2))
    ax.set_xlabel("n. other IS copies in genome")
    sns.despine()


fig, ax = plt.subplots(figsize=(5, 4))
copies_histplot(IS_gains, ax)

plt.ylabel("n. ISs gains")
plt.tight_layout()
plt.savefig(res_fld / "IS_copies.png", dpi=300)
plt.show()
# %%

fig, ax = plt.subplots(figsize=(5, 4))
copies_histplot(is_all, ax)

plt.ylabel("n. all ISs annotations")
plt.tight_layout()
plt.savefig(res_fld / "IS_copies_all.png", dpi=300)
plt.show()

# %%

fig, axs = plt.subplots(1, 2, figsize=(8, 4), sharex=True)

ax = axs[0]
copies_histplot(IS_gains, ax)
ax.set_ylabel("n. ISs gains")
ax.set_title("gained ISs")

ax = axs[1]
copies_histplot(is_all, ax)
ax.set_ylabel("n. all ISs annotations")
ax.set_title("all ISs")

plt.tight_layout()
plt.savefig(res_fld / "IS_copies_combined.png", dpi=300)
plt.show()

# %%
