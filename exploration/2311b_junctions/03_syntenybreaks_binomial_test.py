# %%

import pathlib
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.stats import binom_test


dset = "ST131_ABC"


fld = pathlib.Path(f"../../results/{dset}")
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"

tree = Phylo.read(tree_file, "newick")
tree.ladderize()

L_ter = sum([b.branch_length for b in tree.get_terminals()])
L_int = sum([b.branch_length for b in tree.get_nonterminals()])


print(f"Terminal: {L_ter}")
print(f"Internal: {L_int}")

p = L_ter / (L_ter + L_int)
print(f"Ratio Terminal/Tot: {p}")

# %% events

df_file = f"data/{dset}/singleton_event_df.csv"
dfs = pd.read_csv(df_file, index_col=0)

df_file = f"data/{dset}/mugration/internal_branch_events.csv"
dfi = pd.read_csv(df_file, index_col=0)

df = {
    "internal": {"gain": dfi["gain"].sum(), "loss": dfi["loss"].sum()},
    "terminal": {"gain": dfs["gain"].sum(), "loss": dfs["loss"].sum()},
}
df = pd.DataFrame(df).T
df["tot"] = df["gain"] + df["loss"]
df["corrected loss"] = df["loss"]
df.loc["internal", "corrected loss"] -= 19
df["corrected total"] = df["corrected loss"] + df["gain"]
df.loc["pval"] = np.nan

# %% binomial test gain/losses

for col in df.columns:
    N_int = df[col]["internal"]
    N_ter = df[col]["terminal"]
    pval = binom_test(N_ter, N_int + N_ter, p)
    df.loc["pval", col] = pval

df


# %% synteny breaks
N_in = 4
N_out = 18

# binomial test
p = L_ter / (L_ter + L_int)
print(f"p: {p}")


pval = binom_test(N_out, N_in + N_out, p)
print(f"pval: {pval}")

# %%

# test gain
