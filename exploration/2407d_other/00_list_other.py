# %%
import pandas as pd

fname = "../../results/ST131_ABC/rates/asm20-100-5/terminal_coldspot.csv"
cdf = pd.read_csv(fname)
cdf
# %%

mask = cdf["event_type"] == "other"
cdf = cdf[mask]
cdf
# %%
for e in cdf["edge"]:
    print(e)

# %%
