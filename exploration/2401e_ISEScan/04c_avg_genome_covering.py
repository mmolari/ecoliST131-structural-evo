# %%
import pandas as pd
import numpy as np
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo, SeqIO

fig_fld = pathlib.Path("figs/f04")
fig_fld.mkdir(exist_ok=True, parents=True)

data_fld = pathlib.Path("data/f04")
data_fld.mkdir(exist_ok=True, parents=True)

dset = "ST131_ABC"
fld = pathlib.Path(f"../../results/{dset}")

# genome lengths
iso_L_file = fld / "pangraph/genome_lengths.csv"
iso_L = pd.read_csv(iso_L_file, index_col=0)["length"].to_dict()
# %%


def extract_gbk_product_coverage(gbk_file, pos):
    record = SeqIO.read(gbk_file, "genbank")
    L = len(record.seq)
    cov = np.zeros(L, dtype=bool)
    for feature in record.features:
        if feature.type != "CDS":
            continue
        s, e = int(feature.location.start), int(feature.location.end)
        if s > e:
            s, e = e, s
        if L < 1e5:
            cov[s - 1 : e - 1] = True
        else:
            cov[list(feature.location)] = True
    return cov


df = []
for iso in iso_L.keys():
    print(f"Processing {iso}")
    gbk_file = f"../../data/gbk/{iso}.gbk"
    cov = extract_gbk_product_coverage(gbk_file, iso_L[iso])
    df.append(
        {
            "iso": iso,
            "L": iso_L[iso],
            "product": cov.sum(),
            "empty": (~cov).sum(),
            "mean_cov": cov.mean(),
        }
    )
df = pd.DataFrame(df)
df.to_csv(data_fld / "genome_coverage.csv")

# %%
plt.hist(df["mean_cov"], bins=25)
plt.xlabel("mean genome coverage")
plt.ylabel("n. isolates")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "genome_coverage_hist.png")
plt.show()
# %%
