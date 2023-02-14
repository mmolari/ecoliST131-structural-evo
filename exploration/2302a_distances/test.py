# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib

root_fld = "../.."

fig_fld = pathlib.Path("figs")

df = pd.read_csv(f"{root_fld}/results/ST131/distances/summary-asm20-100-5.csv")
df = df.set_index(["si", "sj"])
mask = df.index.get_level_values(0) > df.index.get_level_values(1)
df = df[mask]

# %%

close_mask = df["core_div"] < 0.0001

# %%
sns.pairplot(df, kind="hist")
plt.savefig(fig_fld / "pairplot.png")
plt.show()

# %%
sns.pairplot(df[close_mask], kind="hist")
plt.savefig(fig_fld / "pairplot_similar.png")
plt.show()
# %%


for yvar in ["mash_dist", "private seq. (bp)", "n. breakpoints", "part. entropy"]:

    sns.histplot(df, x=yvar)
    plt.show()

    sns.histplot(df, y=yvar, x="core_div", cbar=True)
    plt.show()

    sns.histplot(df[close_mask], y=yvar, x="core_div", cbar=True)
    plt.show()

# %%
