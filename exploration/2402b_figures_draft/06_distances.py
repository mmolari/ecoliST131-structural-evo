# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
import pathlib

fig_fld = pathlib.Path("figs/f06")
fig_fld.mkdir(exist_ok=True, parents=True)

fname = "../../results/ST131_ABC/distances/summary-asm20-100-5.csv"
df = pd.read_csv(fname)
mask = df["si"] > df["sj"]
df = df[mask]


with open(
    "../../results/ST131_ABC/pangraph/asm20-100-5-alignment/filtered_corealignment_info_size.json"
) as f:
    aln_info = json.load(f)
core_aln_size = aln_info["core aln size"]
red_aln_size = aln_info["polished aln size"]

df["filt_aln_SNPs"] = df["core_div_filtered"] * red_aln_size

cols = [
    "si",
    "sj",
    "core_div_naive",
    "mash_dist",
    "private seq. (bp)",
    "shared seq. (bp)",
    "n. breakpoints",
    "part. entropy",
    "n. blocks",
    "core_div_filtered",
    "edge_PA",
    "edge_PA_reduced",
    "edge_sharing",
    "block_PA",
    "block_sharing",
    "acc_block_PA",
]


# %%
fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

ax = axs[0]
sns.histplot(data=df, x="filt_aln_SNPs", y="private seq. (bp)", ax=ax)
ax.set_xlabel("SNPs on filtered alignment")

ax = axs[1]
sns.histplot(data=df, x="filt_aln_SNPs", y="block_PA", ax=ax)
ax.set_xlabel("SNPs on filtered alignment")
ax.set_ylabel("block presence/absence distance")


sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "dist_snps.pdf")
plt.show()

# %%
fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))

ax = axs[0]
sns.histplot(data=df, x="core_div_filtered", y="private seq. (bp)", ax=ax)
ax.set_xlabel("divergence on filtered alignment")

ax = axs[1]
sns.histplot(data=df, x="core_div_filtered", y="block_PA", ax=ax)
ax.set_xlabel("divergence on filtered alignment")
ax.set_ylabel("block presence/absence distance")


sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "dist_div.pdf")
plt.show()

# %%
