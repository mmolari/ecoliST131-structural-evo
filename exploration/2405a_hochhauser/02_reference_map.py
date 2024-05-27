# %%
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
from collections import defaultdict

fig_fld = pathlib.Path("figs/n02")
fig_fld.mkdir(exist_ok=True, parents=True)


# read paf file
import pandas as pd


def load_paf(paf_file_path):
    # Define the standard column names according to the PAF format specification
    column_names = [
        "query_name",
        "query_length",
        "query_start",
        "query_end",
        "strand",
        "target_name",
        "target_length",
        "target_start",
        "target_end",
        "residue_matches",
        "alignment_block_length",
        "mapping_quality",
    ]

    try:
        # Read the file into a DataFrame without specifying column names
        df = pd.read_csv(paf_file_path, sep="\t", header=None, comment="#")

        # Separate the standard and optional columns
        standard_df = df.iloc[:, :12]
        optional_df = df.iloc[:, 12:]

        # Assign names to the standard columns
        standard_df.columns = column_names

        # Process the optional fields
        for index, row in optional_df.iterrows():
            for item in row.dropna():
                if ":" in item:
                    key, value = item.split(":", 1)
                    key = key.strip()
                    value = value.strip()
                    # Add the optional field to the standard_df
                    if key not in standard_df.columns:
                        standard_df[key] = pd.NA
                    standard_df.at[index, key] = value

        return standard_df
    except Exception as e:
        print(f"An error occurred while reading the PAF file {fname}:\n{e}")
        return None


fname = "res/map_on_NZ_CP096110.1.paf"
df = load_paf(fname)
df
# %%
fname = "data/hochhauser.csv"
hdf = pd.read_csv(fname, index_col=0)
upstr = hdf.iloc[:, 0].to_dict()
downstr = hdf.iloc[:, 3].to_dict()
hdf = pd.concat([pd.Series(upstr), pd.Series(downstr)], axis=1)
hdf.columns = ["upstream", "downstream"]
hdf
# %%

results = []
for i, row in hdf.iterrows():
    for j in ["upstream", "downstream"]:
        gid = row[j]
        if gid in df["query_name"].values:
            mask = df["query_name"] == gid
            for k, v in df[mask].iterrows():
                if v["residue_matches"] < 0.9 * v["query_length"]:
                    print("warning: short match")
                    continue
                rs = v["target_start"]
                re = v["target_end"]
                results.append(
                    {
                        "hotspot": i,
                        "coregene": j,
                        "gene_id": gid,
                        "start": rs,
                        "end": re,
                    }
                )

results = pd.DataFrame(results)
results.to_csv("res/flanking_genes_loc.csv", index=False)

# %%
lengths = {}
for k, gr in results.groupby("hotspot"):
    if len(gr) == 1:
        continue
    elif len(gr) > 2:
        print("warning: more than 2 flanking genes")
        print(gr)

    assert gr.iloc[0]["end"] > gr.iloc[0]["start"]
    assert gr.iloc[1]["end"] > gr.iloc[1]["start"]
    assert gr.iloc[0]["end"] <= gr.iloc[1]["start"] + 2
    S = min(gr.iloc[0]["end"], gr.iloc[1]["end"])
    E = max(gr.iloc[0]["start"], gr.iloc[1]["start"])
    lengths[k] = E - S + 1

# %%

fig, axs = plt.subplots(
    1, 2, figsize=(9, 5), sharey=True, gridspec_kw={"width_ratios": [5, 1]}
)

ax = axs[0]
for i, row in results.iterrows():
    rs = row["start"]
    re = row["end"]
    u = row["coregene"] == "upstream"
    hs = row["hotspot"]
    marker = "+" if u else "x"
    color = "r" if u else "b"
    ax.plot([rs, re], [hs, hs], color=color, marker=marker)
ax.set_yticks(range(0, len(hdf) + 1, 5))
ax.set_xlabel("NZ_CP096110.1 genome (bp)")
ax.set_ylabel("hotspot number")
ax.grid(axis="y", alpha=0.3)

# legend
ax.plot([], [], marker="+", color="r", ls="", label="upstream")
ax.plot([], [], marker="x", color="b", ls="", label="downstream")
ax.legend()

ax = axs[1]
for k, L in lengths.items():
    ax.barh(k, L, color="k", alpha=0.3)
    ax.plot([1200], [k], color="k", marker=5)
ax.set_xlim(left=1000)
ax.set_xscale("log")
ax.set_xlabel(r"$\Delta L$ (bp)")

sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "flanking_genes_loc.png", dpi=200)
plt.show()

# %%
