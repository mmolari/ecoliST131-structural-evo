# %%
import pandas as pd
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
import utils as ut
from collections import defaultdict

fig_fld = pathlib.Path("figs/f01")
fig_fld.mkdir(parents=True, exist_ok=True)

gene_cluster_json = "../../data/panX/data/ST131_ABC/vis/geneCluster.json"
gc = ut.load_genecluster_json(gene_cluster_json)
gids = gc.index.to_numpy()

loc_fld_real = pathlib.Path("../../results/ST131_ABC/panx/gc_j_pos/asm20-100-5/real")
loc_fld_rand = pathlib.Path("../../results/ST131_ABC/panx/gc_j_pos/asm20-100-5/rand")


def j_per_gene(gids, path):
    ngpg = defaultdict(int)
    for gid in gids:

        row = gc.loc[gid]
        n_loci = len(row["locus"].split())

        loc_file = path / f"genecl_{gid}.csv"
        if not loc_file.exists():
            raise ValueError

        try:
            df = pd.read_csv(loc_file)
            for i, ct in df["id"].value_counts().items():
                ngpg[ct] += 1
            tot = len(df)
            ngpg[0] += n_loci - tot
        except pd.errors.EmptyDataError:
            ngpg[0] += n_loci

    return ngpg


# %%

jpg = j_per_gene(gids, loc_fld_real)
jpg_rand = j_per_gene(gids, loc_fld_rand)

# %%
M = max(max(jpg.keys()), max(jpg_rand.keys()))
x = np.array(range(M + 1))
y = [jpg[i] for i in x]
y_rand = [jpg_rand[i] for i in x]

fig, ax = plt.subplots(1, 1, figsize=(6, 4))
ax.bar(x - 0.2, y, label="real", alpha=0.8, color="C0", width=0.4)
ax.bar(x + 0.2, y_rand, label="random", alpha=0.8, color="C1", width=0.4)
# add numbers on top of bars
ax.bar_label(
    ax.containers[0],
    label_type="edge",
    fontsize=8,
    padding=3,
    fmt="%d",
    rotation=90,
    color="C0",
)
ax.bar_label(
    ax.containers[1],
    label_type="edge",
    fontsize=8,
    padding=3,
    fmt="%d",
    rotation=90,
    color="C1",
)
ax.set_xlabel("n. of junctions associated to the gene")
ax.set_ylabel("n. of genes")
sns.despine()
plt.legend()
plt.tight_layout()
plt.savefig(fig_fld / "j_per_gene.png")
plt.show()

# %%
