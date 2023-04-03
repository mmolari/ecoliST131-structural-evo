# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import numpy as np
import json

fig_fld = pathlib.Path("figs")

df_dupl = pd.read_csv(f"../../results/ST131/distances/summary-asm20-100-5.csv")
df_dupl = df_dupl.set_index(["si", "sj"])
mask = df_dupl.index.get_level_values(0) > df_dupl.index.get_level_values(1)
df_dupl = df_dupl[mask]

df_nodupl = pd.read_csv("../../results/ST131/distances/pangraphasm20-100-5-nodupl.csv")
df_nodupl = df_nodupl.set_index(["si", "sj"])
mask = df_nodupl.index.get_level_values(0) > df_nodupl.index.get_level_values(1)
df_nodupl = df_nodupl[mask]

df = pd.concat([df_nodupl, df_dupl], axis=1, keys=["nodupl", "dupl"])

info_new = "../../results/ST131/pangraph/asm20-100-5-alignment/filtered_corealignment_info.json"
with open(info_new, "r") as f:
    info = json.load(f)
factor = info["polished aln consensus"] + info["polished aln snps"]

df["dupl", "corealn snps"] = df["dupl", "core_div_filtered"] * factor
df["nodupl", "corealn snps"] = df["dupl", "core_div_filtered"] * factor
df["dupl", "private seq. (kbp)"] = df["dupl", "private seq. (bp)"] / 1000
df["dupl", "shared seq. (kbp)"] = df["dupl", "shared seq. (bp)"] / 2000
df["nodupl", "private seq. (kbp)"] = df["nodupl", "private seq. (bp)"] / 1000
df["nodupl", "shared seq. (kbp)"] = df["nodupl", "shared seq. (bp)"] / 2000

# %%
fig, axs = plt.subplots(2, 2, figsize=(8, 8))

ax = axs[0, 0]
xlab, ylab = "corealn snps", "private seq. (kbp)"
sns.histplot(data=df["dupl"], x=xlab, y=ylab, ax=ax)

ax = axs[0, 1]
xlab, ylab = "n. blocks", "private seq. (kbp)"
sns.histplot(data=df["dupl"], x=xlab, y=ylab, ax=ax)

ax = axs[1, 0]
xlab, ylab = "corealn snps", "shared seq. (kbp)"
sns.histplot(data=df["dupl"], x=xlab, y=ylab, ax=ax)

ax = axs[1, 1]
xlab, ylab = "n. blocks", "shared seq. (kbp)"
sns.histplot(data=df["dupl"], x=xlab, y=ylab, ax=ax)

sns.despine(fig)
plt.tight_layout()
# plt.savefig(fig_fld / "pairplot_subset.pdf")
plt.show()
# %%

fig, axs = plt.subplots(2, 4, figsize=(12, 8))

ax = axs[0, 0]
xlab, ylab = "corealn snps", "private seq. (kbp)"
sns.histplot(data=df["nodupl"], x=xlab, y=ylab, ax=ax)

ax = axs[0, 1]
xlab, ylab = "n. blocks", "private seq. (kbp)"
sns.histplot(data=df["nodupl"], x=xlab, y=ylab, ax=ax)

ax = axs[0, 2]
xlab, ylab = "corealn snps", "n. breakpoints"
sns.histplot(data=df["nodupl"], x=xlab, y=ylab, ax=ax)

ax = axs[0, 3]
xlab, ylab = "corealn snps", "part. entropy"
sns.histplot(data=df["nodupl"], x=xlab, y=ylab, ax=ax)

ax = axs[1, 0]
xlab, ylab = "corealn snps", "private seq. (kbp)"
sns.histplot(data=df["dupl"], x=xlab, y=ylab, ax=ax)

ax = axs[1, 1]
xlab, ylab = "n. blocks", "private seq. (kbp)"
sns.histplot(data=df["dupl"], x=xlab, y=ylab, ax=ax)

ax = axs[1, 2]
xlab, ylab = "corealn snps", "n. breakpoints"
sns.histplot(data=df["dupl"], x=xlab, y=ylab, ax=ax)

ax = axs[1, 3]
xlab, ylab = "corealn snps", "part. entropy"
sns.histplot(data=df["dupl"], x=xlab, y=ylab, ax=ax)

sns.despine(fig)
plt.tight_layout()
# plt.savefig(fig_fld / "pairplot_subset.pdf")
plt.show()
# %%

sdf = df["nodupl"]
x = sdf["corealn snps"]
y = sdf["part. entropy"] * sdf["private seq. (kbp)"]
# y = sdf["n. breakpoints"]
# y = sdf["private seq. (kbp)"]
sns.histplot(x=x, y=y)
plt.show()
# %%
