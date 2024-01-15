# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pathlib

# %%

fig_fld = pathlib.Path("figs/f02")
fig_fld.mkdir(exist_ok=True, parents=True)

fname = "../../results/ST131_ABC/annotations/genomad/prophage_summary.tsv"
df = pd.read_csv(fname, sep="\t")


def preprocess(df):
    df["iso"] = df["seq_name"].str.split("|").str[0]
    df["start"] = df["coordinates"].str.split("-").str[0].astype(int)
    df["end"] = df["coordinates"].str.split("-").str[1].astype(int)
    df = df.drop(columns=["fdr", "taxonomy", "topology", "genetic_code", "coordinates"])
    return df


df = preprocess(df)

iso = "NZ_CP124410.1"
# iso = "NZ_CP124326.1"

gnmd_ann_file = f"../../data/genomad/{iso}/{iso}_annotate/{iso}_genes.tsv"
gnmd_ann = pd.read_csv(gnmd_ann_file, sep="\t")
# %%

df_iso = df[df["iso"] == iso]
df_iso

# %%
delta = 8000
proph = df_iso.loc[669]
pp_start = proph["start"]
pp_end = proph["end"]
sqn = proph["seq_name"]
sqn = sqn.replace("provirus_", "")
sqn = sqn.replace("|", "__")

# plot

yt = 1.0
fig, ax = plt.subplots(figsize=(10, 10))

# mark prophage region
ax.axvspan(
    pp_start, pp_end, alpha=0.2, color="royalblue", label="geNomad prophage region"
)
ax.plot([], [], color="gray", marker="|", linewidth=2, label="gene")
ax.plot([], [], color="black", marker="|", linewidth=2, label="prophage hallmark gene")

# add annotations
for i, ann in gnmd_ann.iterrows():
    b, e = ann["start"], ann["end"]
    start_in = (b >= pp_start - delta) and (b <= pp_end + delta)
    end_in = (e >= pp_start - delta) and (e <= pp_end + delta)
    if not (start_in or end_in):
        continue

    color = "black" if ann["virus_hallmark"] else "gray"

    y = 0.5 if ann["strand"] == 1 else 0
    # add wiggle
    y += np.random.uniform(0, 0.5)
    ax.plot([b, e], [y, y], color=color, linewidth=2, marker="|")

    ann_txt = ann["annotation_description"]
    if pd.isna(ann_txt):
        continue
    xt = b
    yt += 0.1
    ax.text(
        xt,
        yt,
        ann["annotation_description"],
        horizontalalignment="left",
        verticalalignment="center",
        fontsize=8,
        color=color,
    )
    ax.plot([b, b], [y, yt], color="gray", linewidth=0.5, linestyle=":")

ax.axhline(0.5, color="gray", linewidth=0.5, linestyle="-")
ax.axhline(1, color="gray", linewidth=0.5, linestyle="-")
ax.set_xlim(pp_start - delta, pp_end + delta)
ax.set_yticks([0.25, 0.75])
ax.set_yticklabels(["-", "+"])
ax.set_ylim(0, yt + 0.1)
ax.set_title(sqn)
ax.set_xlabel("Genome position")
ax.legend(loc="upper left")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / f"{sqn}_annotation.png")
plt.show()

# %%
