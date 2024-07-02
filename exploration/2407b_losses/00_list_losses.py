# %%
import pandas as pd

fname = "../../results/ST131_ABC/rates/asm20-100-5/terminal_coldspot.csv"
cdf = pd.read_csv(fname)
cdf
# %%

mask = cdf["event_type"] == "loss"
cdf[mask]
# %%
for e in cdf[mask]["edge"].to_list():
    print(e)

# %%
