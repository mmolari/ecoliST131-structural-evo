# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
import scipy.stats as stats

fig_fld = pathlib.Path("figs/f11")
fig_fld.mkdir(parents=True, exist_ok=True)

Ls = "../../results/ST131_ABC/pangraph/genome_lengths.csv"
Ls = pd.read_csv(Ls, index_col=0)
metadata = "../../results/ST131_ABC/metadata.csv"
metadata = pd.read_csv(metadata, index_col=0)
Ls["year"] = metadata["year"]

# %%
R = sns.regplot(
    data=Ls,
    x="year",
    y="length",
    scatter_kws={"alpha": 0.5},
)
# linear regression
mask = ~Ls.isna().any(axis=1)
slope, intercept, r_value, p_value, std_err = stats.linregress(
    Ls[mask]["year"], Ls[mask]["length"]
)
# display line equation
R.text(0.1, 0.9, f"y = {slope:.0f} x + {intercept:.2e}", transform=R.transAxes)

plt.xlabel("Year")
plt.ylabel("Genome length (bp)")
plt.savefig(fig_fld / "genome_length_vs_year.png")
plt.show()


# %%
