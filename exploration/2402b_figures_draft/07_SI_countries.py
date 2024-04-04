# %%

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pathlib

# assign a color to each continent
continent_color = {
    "Africa": sns.color_palette("Set2")[3],
    "Asia": sns.color_palette("Set2")[2],
    "Europe": sns.color_palette("Set2")[1],
    "North America": sns.color_palette("Set2")[0],
    "Oceania": sns.color_palette("Set2")[5],
    "South America": sns.color_palette("Set2")[4],
}


# %%

fname = "../../results/ST131_ABC/metadata.csv"
df_met = pd.read_csv(fname, index_col=0)
svfld = pathlib.Path("figs/f07")
svfld.mkdir(parents=True, exist_ok=True)

# %%
fig, axs = plt.subplots(2, 1, figsize=(8, 7))

ax = axs[0]
sns.countplot(
    data=df_met,
    x="geo_loc_name",
    ax=ax,
    hue="continent",
    palette=continent_color,
    order=df_met["geo_loc_name"].value_counts().index,
)
ax.set_xticks(range(len(df_met["geo_loc_name"].value_counts())))
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
cts = df_met["continent"].value_counts().to_dict()
ax.set_ylabel("number of isolates")
ax.set_xlabel("")
ax.set_title("isolates per country")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=[f"{l} (n={cts[l]})" for l in labels])
ax.grid(axis="y", linestyle="--", alpha=0.5)

df_met["year_cat"] = pd.cut(
    df_met["year"],
    [2000, 2010, 2020, 2030],
    labels=["2000-2009", "2010-2019", "2020-2022"],
    right=False,
)

year_colors = {
    "2000-2009": sns.color_palette("Blues")[1],
    "2010-2019": sns.color_palette("Blues")[3],
    "2020-2022": sns.color_palette("Blues")[5],
}

ax = axs[1]
sns.histplot(
    data=df_met,
    x="year",
    ax=ax,
    discrete=True,
    color="lightgrey",
    edgecolor="black",
    hue="year_cat",
    palette=year_colors,
)
# redefine legend
handles = [
    mpl.patches.Patch(facecolor=year_colors[k], alpha=0.5, edgecolor="k")
    for k in year_colors
]
labels = [
    f"{k} (n={v})" for k, v in df_met["year_cat"].value_counts().sort_index().items()
]
ax.legend(handles=handles, labels=labels)
ax.set_ylabel("number of isolates")
ax.set_xlabel("")
ax.set_title("isolates per year")
ax.grid(axis="y", linestyle="--", alpha=0.5)

sns.despine()
plt.tight_layout()
fig.savefig(svfld / "metadata_hist.pdf")
plt.show()
# %%
