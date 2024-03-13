# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib

fig_fld = pathlib.Path("figs/n4")
fig_fld.mkdir(exist_ok=True, parents=True)


def load_df():
    fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
    df = pd.read_csv(fname, index_col=0)

    fnames = {
        "df": "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/defensefinder_real.csv",
        "gm": "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/genomad_real.csv",
        "if": "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/integronfinder_real.csv",
        "is": "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/ISEScan_real.csv",
    }

    for k, fname in fnames.items():
        df2 = pd.read_csv(fname, index_col=0)
        df2 = df2["junction"].value_counts()
        df[f"{k}"] = df.index.map(df2)
        df[f"{k}"] = df[f"{k}"].fillna(0)
        df[f"{k}"] = df[f"{k}"].astype(int)
        df[f"{k}"] = df[f"{k}"] > 0

    return df


js = load_df()
js
# %%
mask = js["n_categories"] == 2
print(f"n. two categories = {mask.sum()} / {mask.size}")
js = js[mask]
# %%
assert np.all(js["singleton"] == (js["n_iso"] - 1 == js["majority_category"]))
js["n_iso"].value_counts()
# %%
mask = js["n_iso"] == 222
js = js[mask]
js["singleton"].value_counts()
# %%
# add MGE annotations
# add internal vs terminal
js[["gm", "if", "is"]].value_counts()
js["cat"] = "other"

for k, val in [
    ("IS", (False, False, True)),
    ("none", (False, False, False)),
    ("prophage", (True, False, True)),
    ("prophage", (True, False, False)),
    ("integron", (False, True, True)),
]:
    keys = ["gm", "if", "is"]
    mask = np.all(js[keys] == val, axis=1)
    js.loc[mask, "cat"] = k
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
ax.set_title("backbone coldspots breakdonw")

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / f"backbone_coldspots_breakdonw.png")
plt.show()

# %%
