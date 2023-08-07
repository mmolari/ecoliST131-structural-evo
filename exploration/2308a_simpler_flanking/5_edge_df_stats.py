# %%
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

import utils as ut

N_strains = len(ut.load_tree().get_terminals())

fig_fld = ut.fig_fld / "core_edges"
fig_fld.mkdir(exist_ok=True)
# %%


df = pd.read_csv(ut.expl_fld / "core_junctions_df.csv", index_col=0)
df.index = [ut.Edge.from_str_id(e) for e in df.index]
df["delta_len"] = df["len_max"] - df["len_min"]
# %%
fig, axs = plt.subplots(1, 3, figsize=(12, 4))
ax = axs[0]
sns.scatterplot(
    data=df,
    x="avg_len",
    y="avg_n_blocks",
    hue="n_strains",
    hue_norm=mpl.colors.Normalize(vmin=1, vmax=N_strains),
    palette="rainbow",
    alpha=0.5,
    ax=ax,
)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("avg. junction length (bp)")
ax.set_ylabel("avg. n. blocks in junction")
ax.grid(True, alpha=0.3)

ax = axs[1]
sns.scatterplot(
    data=df,
    x="len_max",
    y="len_min",
    hue="n_strains",
    hue_norm=mpl.colors.Normalize(vmin=1, vmax=N_strains),
    palette="rainbow",
    alpha=0.5,
    ax=ax,
)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("max junction len (bp)")
ax.set_ylabel("min junction length (bp)")
ax.grid(True, alpha=0.3)

ax = axs[2]
sns.scatterplot(
    data=df,
    x="avg_len",
    y="n_strains",
    hue="avg_n_blocks",
    hue_norm=mpl.colors.LogNorm(),
    palette="rainbow",
    alpha=0.5,
    ax=ax,
)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("avg. junction length (bp)")
ax.set_ylabel("n. strains")
ax.grid(True, alpha=0.3)


sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "edge_stats.png", facecolor="w")
plt.show()

# %%
