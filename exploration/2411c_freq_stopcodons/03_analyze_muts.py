# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pathlib
import utils as ut

ld1 = pathlib.Path("res/f00")
ld2 = pathlib.Path("res/f02")
fig_fld = pathlib.Path("res/f03")
fig_fld.mkdir(exist_ok=True)


# load gene clusters info
gene_cluster_json = "../../data/panX/data/ST131_ABC/vis/geneCluster.json"
gdf = ut.load_genecluster_json(gene_cluster_json)


def load_files(fld):
    Mc = pd.read_csv(fld / "mut_count.csv", index_col=[0, 1, 2])  # mutation counts
    Pm = pd.read_csv(fld / "mut_possible.csv", index_col=[0, 1, 2])  # total counts
    As = pd.read_csv(fld / "aln_stats.csv", index_col=0).T  # background counts

    return {
        "mut_count": Mc,
        "mut_possible": Pm,
        "aln_stats": As,
    }


# %%
core_files = load_files(ld1)
almost_files = load_files(ld2)

# %%

MCc = core_files["mut_count"]
MCa = almost_files["mut_count"]

MC = MCc.join(MCa, how="outer").fillna(0).astype(int)
MC = MC.sum(axis=1).sort_values()
MC = MC.reset_index().rename(
    columns={"level_0": "from", "level_1": "to", "level_2": "type", 0: "count"}
)
tot_muts_mc = MC.groupby(["from", "to"])["count"].sum()
# %%

PMc = core_files["mut_possible"]
PMa = almost_files["mut_possible"]

PM = PMc.join(PMa, how="outer").fillna(0).astype(int)
PM = PM.sum(axis=1).sort_values()
PM = PM.reset_index().rename(
    columns={"level_0": "from", "level_1": "to", "level_2": "type", 0: "count"}
)
PM = PM.pivot_table(
    index=["from", "to"], columns="type", values="count", fill_value=0
).astype(int)


# %%
Pm_frac = PM / PM.sum(axis=1).values[:, None]
exp_stops = Pm_frac["*"] * tot_muts_mc
obs_stops = MC[MC["type"] == "*"].set_index(["from", "to"])["count"]

stops = pd.DataFrame({"obs": obs_stops, "exp": exp_stops}).fillna(0)
# filter out rows that are not completely zero
stops = stops[(stops.T != 0).any()]
stops.reset_index(inplace=True)
stops["mutation"] = stops["from"] + "->" + stops["to"]
stops["frac"] = stops["obs"] / stops["exp"]
stops
# %%
# plot
fig, ax = plt.subplots(figsize=(6, 5))
sns.barplot(data=stops, x="mutation", y="frac", ax=ax)
# annotated number of observed over expected
for i, row in stops.iterrows():
    ax.text(
        i, row["frac"], f"{row['obs']:.0f}/{row['exp']:.0f}", ha="center", va="bottom"
    )
tot_obs = stops["obs"].sum()
tot_exp = stops["exp"].sum()

plt.axhline(
    tot_obs / tot_exp,
    ls="--",
    color="silver",
    label=f"average: {tot_obs:.0f}/{tot_exp:.0f}",
    zorder=-2,
)
plt.legend()
# plt.xticks(rotation=90)
plt.ylabel("fraction of observed/expected nonsense mutations")
sns.despine()
plt.tight_layout()
plt.savefig("res/f03/obs_exp_nonsense_mutations.png")
plt.show()

# %%
