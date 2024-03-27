# %%
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gainloss_utils as ut
import pathlib

fig_fld = pathlib.Path("figs/f05")
fig_fld.mkdir(exist_ok=True, parents=True)


def svfig(svname):
    pass
    plt.savefig(fig_fld / f"{svname}.png", dpi=200)
    plt.savefig(fig_fld / f"{svname}.svg")


# %% breakdown

js = ut.load_cdf(ut.fnames["coldspots"])
fig, axs = plt.subplots(
    1,
    2,
    figsize=(7, 4),
    gridspec_kw={"width_ratios": [2, 1]},
    sharey=True,
)

ax = axs[0]
et = js["cat"].value_counts().sort_index()
x = 0
xticks, xlabels = [], []
for k, v in et.items():
    b = ax.bar(x, v, color=ut.cat_colors[k])
    ax.bar_label(b)
    xticks.append(x)
    xlabels.append(k)
    x += 1

ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)
ax.set_ylabel("n. coldspots")

# singletons
ax = axs[1]
colors = {
    True: "gray",
    False: "lightgray",
}

et = js["singleton"].value_counts().to_dict()
xticks, xlabels = [], []
x = 0
for k, v in et.items():
    b = ax.bar(x, v, color=colors[k])
    ax.bar_label(b)
    xticks.append(x)
    xlabels.append("singleton" if k else "non-singleton")
    x += 1
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)
# ax.set_ylabel("n. coldspots")

sns.despine()
plt.tight_layout()
svfig("suppl_backbone_coldspots_breakdonw")
plt.show()

# %%
