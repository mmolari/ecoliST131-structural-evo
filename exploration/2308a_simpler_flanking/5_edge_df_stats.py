# %%
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

import utils as ut

N_strains = len(ut.load_tree().get_terminals())
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
    x="avg_len",
    y="delta_len",
    hue="avg_n_blocks",
    hue_norm=mpl.colors.LogNorm(),
    palette="rainbow",
    alpha=0.5,
    ax=ax,
)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("avg. junction length (bp)")
ax.set_ylabel("max - min junction length (bp)")
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
plt.show()


# %%

sns.scatterplot(
    data=df,
    x="len",
    y="len_var",
    # hue="len",
    # palette="viridis",
    # hue_norm=mpl.colors.LogNorm(),
    alpha=0.5,
)
plt.xscale("log")
plt.yscale("log")
plt.tight_layout()
plt.show()

sns.scatterplot(
    data=df,
    x="len",
    y="len_dupl",
    # hue="len",
    # palette="viridis",
    # hue_norm=mpl.colors.LogNorm(),
    alpha=0.5,
)
plt.xscale("log")
plt.yscale("log")
plt.tight_layout()
plt.show()

sns.scatterplot(
    data=df,
    x="len",
    y="n_strains",
    # hue="len",
    # palette="viridis",
    # hue_norm=mpl.colors.LogNorm(),
    alpha=0.5,
)
plt.xscale("log")
plt.yscale("log")
plt.tight_layout()
plt.show()
