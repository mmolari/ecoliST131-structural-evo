# %%
import pandas as pd
import seaborn as sns


fname = "../../results/ST131_ABC/distances/summary-asm20-100-5.csv"
df = pd.read_csv(fname)
df

# %%
mask = df["si"] > df["sj"]
df[mask]["private seq. (bp)"].mean()

# %%
