# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import numpy as np
import json

root_fld = "../.."

fig_fld = pathlib.Path("figs/scatter")
fig_fld.mkdir(exist_ok=True)

df = pd.read_csv(f"{root_fld}/results/ST131/distances/summary-asm20-100-5.csv")
df = df.set_index(["si", "sj"])
mask = df.index.get_level_values(0) > df.index.get_level_values(1)
df = df[mask]

info_new = "../../results/ST131/pangraph/asm20-100-5-alignment/filtered_corealignment_info.json"
with open(info_new, "r") as f:
    info = json.load(f)
factor = info["polished aln consensus"] + info["polished aln snps"]

df["corealn snps"] = df["core_div_filtered"] * factor
df["private seq. (kbp)"] = df["private seq. (bp)"] / 1000
df["shared seq. (kbp)"] = df["shared seq. (bp)"] / 2000

# %%
fig, axs = plt.subplots(2, 2, figsize=(8, 8))

ax = axs[0, 0]
xlab, ylab = "corealn snps", "private seq. (kbp)"
sns.histplot(data=df, x=xlab, y=ylab, ax=ax)
# x,y =df[xlab].to_numpy(), df[ylab].to_numpy()
# m, q = np.polyfit(x,y,deg=1)
# a,b = x.max(), x.min()
# ax.plot([a,b], [a*m + q, b*m + q], "k--")

ax = axs[0, 1]
xlab, ylab = "n. blocks", "private seq. (kbp)"
sns.histplot(data=df, x=xlab, y=ylab, ax=ax)
x, y = df[xlab].to_numpy(), df[ylab].to_numpy()
m, q = np.polyfit(x, y, deg=1)
# m, q = np.sum(x*y)/np.sum(x*x), 0
a, b = x.max(), x.min()
ax.plot([a, b], [a * m + q, b * m + q], "k--")
ax.set_title(f"y = {m:.1f} x + {q:.0f}")

ax = axs[1, 0]
xlab, ylab = "corealn snps", "shared seq. (kbp)"
sns.histplot(data=df, x=xlab, y=ylab, ax=ax)
# x,y =df[xlab].to_numpy(), df[ylab].to_numpy()
# m, q = np.polyfit(x,y,deg=1)
# a,b = x.max(), x.min()
# ax.plot([a,b], [a*m + q, b*m + q], "k--")

ax = axs[1, 1]
xlab, ylab = "n. blocks", "shared seq. (kbp)"
sns.histplot(data=df, x=xlab, y=ylab, ax=ax)
# x, y = df[xlab].to_numpy(), df[ylab].to_numpy()
# m, q = np.polyfit(x, y, deg=1)
# a, b = x.max(), x.min()
# ax.plot([a, b], [a * m + q, b * m + q], "k--")
# ax.set_title(f"y = {m:.1f} x + {q:.1e}")

sns.despine(fig)
plt.tight_layout()
plt.savefig(fig_fld / "pairplot_subset.pdf")
plt.show()
# %%

pairs = [
    ("corealn snps", "private seq. (bp)", "privseq"),
    ("corealn snps", "n. blocks", "nblocks"),
    ("corealn snps", "mash_dist", "mash"),
]

for xl, yl, sl in pairs:
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sns.histplot(data=df, x=xl, y=yl, ax=ax)
    sns.despine()
    plt.tight_layout()
    plt.savefig(fig_fld / f"hist_{sl}.png")
    plt.show()

# %%
