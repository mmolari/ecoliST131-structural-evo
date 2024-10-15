# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import pathlib
from collections import defaultdict
import pypangraph as pp


fld = pathlib.Path("../../results/ST131_ABC")
# fig_fld = pathlib.Path("figs/f00")
# fig_fld.mkdir(exist_ok=True, parents=True)

# colors_file = fld / "pangraph/coresynt-asm20-100-5/blocks.csv"
# mergers_file = fld / "pangraph/coresynt-asm20-100-5/mergers.csv"

# cl_df = pd.read_csv(colors_file, index_col=0)
# mg_df = pd.read_csv(mergers_file, index_col=0)

pangraph_file = fld / "pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()


# %%
fig, ax = plt.subplots(figsize=(10, 5))
sns.histplot(
    data=bdf,
    x="n. strains",
    weights="len",
    discrete=True,
    ax=ax,
    cumulative=True,
)
plt.tight_layout()
plt.show()


# %%
bdf[bdf["n. strains"] < 20]["len"].sum() / 1e6
# %%
