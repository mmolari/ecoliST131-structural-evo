# %%
import pandas as pd
import json
from collections import defaultdict


# %%
fname = "../../config/datasets/ST131_full/acc_info.csv"
df = pd.read_csv(fname)

# %%
# extract entry of type "chromosome" for each assembly
is_chr = df["type"] == "chromosome"
chr_ids = df[is_chr].set_index("asmbl")["id"].to_dict()
# %%

pl_ids = defaultdict(list)
for idx, row in df[~is_chr].iterrows():
    chromosome = chr_ids[row["asmbl"]]
    pl_ids[chromosome].append(row["id"])
pl_ids
# %%
out_fname = "data/plasmids.json"
with open(out_fname, "w") as handle:
    json.dump(pl_ids, handle, indent=4)
# %%
