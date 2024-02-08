# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib

fig_fld = pathlib.Path("figs/f01")
fig_fld.mkdir(exist_ok=True, parents=True)

fld = pathlib.Path("../../results/ST131_ABC")

el_fname = fld / "annotations/junct_pos_asm20-100-5/ISEScan_real.csv"
df_el = pd.read_csv(el_fname)

el_rand_fname = fld / "annotations/junct_pos_asm20-100-5/ISEScan_rand.csv"
df_el_rand = pd.read_csv(el_rand_fname)

edge_fname = fld / "backbone_joints/asm20-100-5/edge_len.csv"
df_edge = pd.read_csv(edge_fname, index_col=0)
edges = df_edge.columns.to_list()

loc_fname = fld / "annotations/loc/ISEScan.csv"
df_el_loc = pd.read_csv(loc_fname, index_col=0)
# %%

# IS per junction

IS_per_J = pd.DataFrame(edges, columns=["junctions"]).set_index(
    "junctions", verify_integrity=True
)
IS_per_J["real"] = df_el["junction"].value_counts()
IS_per_J["random"] = df_el_rand["junction"].value_counts()
IS_per_J.fillna(0, inplace=True)

fig, axs = plt.subplots(2, 1, figsize=(8, 7))

for lab, c in [("real", "C0"), ("random", "gray")]:
    for nax, w in enumerate([None, IS_per_J[lab]]):
        ax = axs[nax]
        sns.histplot(
            x=IS_per_J[lab],
            weights=w,
            ax=ax,
            label=lab,
            element="step",
            cumulative=True,
            fill=False,
            stat="probability",
            discrete=True,
            color=c,
        )

ax = axs[0]
ax.set_xlabel("n. of IS in junction")
ax.set_ylabel("fraction of junctions")
ax.set_xscale("symlog", linthresh=1)
ax.set_ylim(0, 1.02)
ax.legend()


ax = axs[1]
ax.set_xlabel("n. of IS in junction")
ax.set_ylabel("fraction of IS")
ax.set_ylim(0, 1.02)
ax.legend()

plt.tight_layout()
plt.savefig(fig_fld / "IS_per_junction.png")
plt.show()


# %%
# junction per IS
J_per_IS = pd.DataFrame(index=df_el_loc.index)
J_per_IS["real"] = df_el.groupby("id")["junction"].nunique()
J_per_IS["random"] = df_el_rand.groupby("id")["junction"].nunique()
J_per_IS.fillna(0, inplace=True)
J_per_IS
# %%


fig, ax = plt.subplots(1, 1, figsize=(8, 4))

sns.histplot(J_per_IS["real"], ax=ax, label="real", discrete=True)
sns.histplot(J_per_IS["random"], ax=ax, color="gray", label="random", discrete=True)
ax.set_xlabel("n. of junctions")
ax.set_ylabel("n. of IS")
ax.legend()

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "junction_per_IS.png")
plt.show()

# %%
