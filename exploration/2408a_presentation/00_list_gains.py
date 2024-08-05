# %%
import pandas as pd

fname = "../../results/ST131_ABC/rates/asm20-100-5/terminal_coldspot.csv"
fname = "../../results/ST131_ABC/rates/asm20-100-5/nonsingleton_junct_df.csv"
cdf = pd.read_csv(fname)
cdf
# %%

# mask = cdf["event_type"] == "gain"
# cdf = cdf[mask]
# cdf
# %%
# lab = "defensefinder"
lab = "isescan"
# lab = "genomad"
mask = cdf[lab] > 0
sdf = cdf[mask].sort_values("majority_category", ascending=False)
sdf

# %%
