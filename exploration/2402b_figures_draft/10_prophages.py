# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib

fig_fld = pathlib.Path("figs/f10")
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
df

# %%

fig, axs = plt.subplots(1, 2, figsize=(7, 3))
ax = axs[0]
sns.histplot(
    df["iso"].value_counts(),
    discrete=True,
    ax=ax,
)
xticks = range(0, df["iso"].value_counts().max() + 1, 2)
ax.set_xticks(xticks)
ax.set_xlabel("n. of prophages per genome")
ax.set_ylabel("n. of isolates")

ax = axs[1]
sns.histplot(df["length"], ax=ax)
ax.set_xlabel("prophage length [bp]")
ax.set_ylabel("n. of prophages")

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "prophage_summary.pdf")
plt.show()

# %%
