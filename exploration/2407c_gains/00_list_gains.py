# %%
import pandas as pd

fname = "../../results/ST131_ABC/rates/asm20-100-5/terminal_coldspot.csv"
cdf = pd.read_csv(fname)
cdf
# %%

mask = cdf["event_type"] == "gain"
cdf = cdf[mask]
cdf
# %%

for lab in ["genomad", "integrons", "isescan", "defensefinder"]:
    mask = cdf[lab] > 0

    # randomly pick two elements
    N = min(3, mask.sum())
    for e in cdf[mask]["edge"].sample(N).to_list():
        print(e)


# %%
