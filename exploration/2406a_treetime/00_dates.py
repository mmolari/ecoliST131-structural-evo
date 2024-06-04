# %%
import pandas as pd

df = pd.read_csv("../../config/datasets/ST131_ABC/metadata.csv", index_col=0)
sdf = df[["year"]].copy()
sdf.index.name = "name"
# format year as YYYY-XX-XX since month and day are unknown
sdf["year"] = sdf.apply(
    lambda x: f"{x['year']:.0f}-XX-XX" if not pd.isna(x["year"]) else None, axis=1
)
sdf.rename(columns={"year": "date"}, inplace=True)
sdf.to_csv("data/dates.csv")
# %%
