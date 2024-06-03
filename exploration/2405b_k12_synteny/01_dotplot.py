# %%
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
from dataclasses import dataclass

fig_fld = pathlib.Path("figs/n01")
fig_fld.mkdir(exist_ok=True, parents=True)

h = "qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen slen"

df = pd.read_csv("data/blast_res/k12_on_ref.tsv", sep="\t", header=None)
df.columns = h.split()
df["div"] = df["mismatch"] / (df["length"] - df["gaps"])
df
# %%

sns.histplot(
    data=df, x="div", weights="length", bins=100, stat="density", element="step"
)
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "div.png")
plt.show()
# %%

cmap = mpl.colormaps.get("coolwarm_r")
norm = mpl.colors.Normalize(vmax=100, vmin=90, clip=True)
fig, ax = plt.subplots(1, 1, figsize=(8, 8))

mask = df["length"] >= 5000

sdf = df[mask]


avg_pident = (sdf["div"] * sdf["length"] / sdf["length"].sum()).sum()
print(f"{avg_pident=}")

for idx, row in sdf.iterrows():
    qs, qe = row[["qstart", "qend"]]
    rs, re = row[["sstart", "send"]]
    # assert qs < qe
    # assert rs < re, f"{rs=}, {re=}"
    # color = cmap(norm(row["pident"]))
    color = "k"
    ax.plot([qs, qe], [rs, re], color=color)

# ax.set_title(f"avg. identity = {avg_pident:.1f} %")
ax.set_xlabel(f"E.coli K12")
ax.set_ylabel(df["sseqid"].unique()[0])
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(fig_fld / "k12_dotplot.png")
plt.savefig(fig_fld / "k12_dotplot.pdf")
plt.show()

# %%


@dataclass
class Interval:
    s: int
    e: int


def rectify_interval(s, e):
    I = Interval(s, e) if s < e else Interval(e, s)
    return I


def has_overlap(I1, I2):
    assert I1.s < I1.e
    assert I2.s < I2.e
    if I1.e < I2.s:
        return False
    elif I2.e < I1.s:
        return False
    return True


mask = df["length"] >= 500
sdf = df[mask]

fig, ax = plt.subplots(1, 1, figsize=(8, 8))
forbidden = set()
N = len(sdf)
for i in range(N):
    row = sdf.iloc[i]
    qs, qe = row["qstart"], row["qend"]
    rs, re = row["sstart"], row["send"]
    Iq = rectify_interval(qs, qe)
    Ir = rectify_interval(rs, re)
    for j in range(i + 1, N):
        row_j = sdf.iloc[j]
        Iqj = rectify_interval(row_j["qstart"], row_j["qend"])
        Irj = rectify_interval(row_j["sstart"], row_j["send"])
        if has_overlap(Iq, Iqj) or has_overlap(Ir, Irj):
            forbidden.update([i, j])
    if i in forbidden:
        continue
    ax.plot([qs, qe], [rs, re], color="k")


ax.set_xlabel(f"E.coli K12")
ax.set_ylabel(df["sseqid"].unique()[0])
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(fig_fld / "k12_dotplot_unique.png")
plt.savefig(fig_fld / "k12_dotplot_unique.pdf")
plt.show()

# %%
