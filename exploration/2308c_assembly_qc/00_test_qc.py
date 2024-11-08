# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

df = pd.read_csv("../../results/ST131/assembly_qc/summary.csv", index_col="iso")

# %%
columns = [
    "one_line_summary",
    "Complete",
    "Single copy",
    "Multi copy",
    "Fragmented",
    "Missing",
    "n_markers",
    "domain",
    "Number of scaffolds",
    "Number of contigs",
    "Total length",
    "Percent gaps",
    "Scaffold N50",
    "Contigs N50",
]
# %%
sns.scatterplot(data=df, x="Contigs N50", y="Percent gaps")

# %%

fig, ax = plt.subplots(figsize=(10, 15))

colors = {
    "Single copy": "green",
    "Multi copy": "darkgreen",
    "Fragmented": "orange",
    "Missing": "red",
}

y = 0
for iso, row in df.iterrows():
    x = 0
    tot = sum(row[k] for k in colors.keys())
    for k, c in colors.items():
        if k in ["Single copy", "Multi copy"]:
            continue
        dx = (row[k] / tot) * 100
        ax.barh(y, dx, left=x, color=c, label=k if y == 0 else None)
        x += dx
    y += 1
ax.legend()
ax.set_yticks(range(len(df)))
ax.set_ylim(-1, len(df))
ax.set_yticklabels(df.index)
# ax.axvline(2)
sns.despine()
ax.set_xlabel("percent of missing core genes")
plt.tight_layout()
# plt.savefig("figs/missing_core_genes.png", facecolor="white")
plt.show()
# %%
fig, ax = plt.subplots(figsize=(5, 4))
sns.histplot(data=df, x="Number of contigs", binwidth=1, ax=ax)
for iso, row in df.iterrows():
    nc = row["Number of contigs"]
    if nc > 1:
        ax.annotate(
            f"{iso} ({nc})",
            xy=(nc, 2),
            xytext=(nc - 10, np.random.randint(10, 30)),
            arrowprops=dict(arrowstyle="->"),
        )
sns.despine()
plt.tight_layout()
# plt.savefig("figs/n_contigs.png", facecolor="white")
plt.show()

# %%
