# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import re

fig_fld = pathlib.Path("figs/f02")
fig_fld.mkdir(parents=True, exist_ok=True)

fnames = [
    "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/genomad_real.csv",
    "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/integronfinder_real.csv",
    "../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/ISEScan_real.csv",
]

junct_file = "../../results/ST131_ABC/backbone_joints/asm20-100-5/edge_count.csv"
Js = pd.read_csv(junct_file, index_col=0)
Js.rename_axis("junction", inplace=True)

junct_info_file = (
    "../../results/ST131_ABC/backbone_joints/asm20-100-5/junctions_stats.csv"
)
Jstats = pd.read_csv(junct_info_file, index_col=0)

dfs = {}
for fname in fnames:
    f = pathlib.Path(fname)
    df = pd.read_csv(f)
    tool = re.search(r"(\w+)_real.csv", f.name).group(1)
    dfs[tool] = df

for name, df in dfs.items():
    Jct = df["junction"].value_counts()
    Js[name] = Js.index.map(Jct).fillna(0).astype(int)


Js = Js.join(Jstats)
Js["delta_len"] = Js["max_length"] - Js["min_length"]


# %%

for tool in dfs:
    g = sns.JointGrid(
        Js,
        x="delta_len",
        y="n_categories",
        hue=Js[tool] > 0,
    )
    g.plot_joint(sns.scatterplot, alpha=0.2)
    g.plot_marginals(
        sns.histplot, log_scale=True, bins=25, multiple="stack", linewidth=0.5
    )
    g.ax_joint.set_xlabel(r"junction $\Delta L$")
    g.ax_joint.set_ylabel("n. of path categories")
    plt.tight_layout()
    plt.savefig(fig_fld / f"joint_{tool}.png")
    plt.show()

# %%
Js["MGEs"] = ""
Js.loc[Js["ISEScan"] > 0, "MGEs"] += "|is"
Js.loc[Js["genomad"] > 0, "MGEs"] += "|pp"
# Js.loc[Js["integronfinder"] > 0, "MGEs"] += "|in"
Js["MGEs"] = Js["MGEs"].str.lstrip("|")
Js["MGEs"].replace("", "none", inplace=True)
Js["MGEs"] = pd.Categorical(
    Js["MGEs"],
    # categories=["none", "is", "pp", "in", "is|pp", "is|in", "pp|in", "is|pp|in"],
    categories=["none", "is", "pp", "is|pp"],
)
Js["MGEs"].value_counts(dropna=True)
# %%
g = sns.JointGrid(
    Js,
    x="delta_len",
    y="n_categories",
    hue="MGEs",
    height=6,
    ratio=4,
)
g.plot_joint(sns.scatterplot, alpha=0.3)
g.plot_marginals(
    sns.histplot,
    log_scale=True,
    bins=25,
    multiple="stack",
    element="step",
    linewidth=0.4,
)
g.ax_joint.set_xlabel(r"junction $\Delta L$")
g.ax_joint.set_ylabel("n. of path categories")
plt.tight_layout()
plt.savefig(fig_fld / "joint_MGEs.png")
plt.show()
# %%
