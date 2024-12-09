# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
import pathlib
import itertools as itt

fig_fld = pathlib.Path("figs/f06")
fig_fld.mkdir(exist_ok=True, parents=True)

fname = "res/blockdist.csv"
bdf = pd.read_csv(fname)
mask = bdf["si"] > bdf["sj"]
bdf = bdf.loc[mask, :]

fname = "res/pwdist.csv"
pdf = pd.read_csv(fname)
mask = pdf["si"] > pdf["sj"]
pdf = pdf.loc[mask, :]

# %%
sns.histplot(data=pdf, x="private seq. (bp)", bins=50)
plt.show()

sns.histplot(data=bdf, x="block_PA", bins=50)
plt.show()


# %%
