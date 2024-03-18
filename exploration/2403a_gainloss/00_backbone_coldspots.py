# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import utils as ut

fig_fld = pathlib.Path("figs/n0")
fig_fld.mkdir(exist_ok=True, parents=True)

fname = "../../results/ST131_ABC/rates/asm20-100-5/coldspot_df.csv"
js = pd.read_csv(fname, index_col=0)
js
# %%
# add MGE annotations
js = ut.assign_category(js)
js["cat"].value_counts()
# %%
js["singleton"].value_counts()
# %%

fig, ax = plt.subplots(1, 1, figsize=(4, 4))

colors = {
    "IS": "C0",
    "prophage": "C4",
    "integron": "C1",
    "none": "#b8b8b8",
}
et = js["cat"].value_counts().to_dict()
tot = sum(et.values())
bottom = 0
for k, v in et.items():
    ax.bar(0, v, 0.8, bottom=bottom, label=k, color=colors[k])
    ax.text(0, bottom + v / 2, f"{k} (n={v})", ha="center", va="center")
    bottom += v

# singletons
bottom = 0
colors = {
    True: "gray",
    False: "lightgray",
}

et = js["singleton"].value_counts().to_dict()
for k, v in et.items():
    ax.bar(1, v, 0.8, bottom=bottom, label=k, color=colors[k])
    lab = "singleton" if k else "non-singleton"
    ax.text(1, bottom + v / 2, f"{lab} (n={v})", ha="center", va="center")
    bottom += v

ax.set_xticks([0, 1])
ax.set_xticklabels(["annotations", "frequency"])
ax.set_title("backbone coldspots breakdown")

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / f"backbone_coldspots_breakdonw.png")
plt.show()

# %%
