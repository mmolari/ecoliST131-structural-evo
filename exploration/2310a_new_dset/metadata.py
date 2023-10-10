# %%
import pandas as pd
import pycountry


md_file = "../../config/new_dset/n225-metadata-cleaned.csv"
df = pd.read_csv(md_file)


def normalize_country_name(name):
    try:
        return pycountry.countries.lookup(name).name
    except:
        print("not matched: ", name)
        return None


df["location"] = df["country"].apply(normalize_country_name)

# %%
for l in df["location"]:
    print(l)
# %%

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
country_counts = df["location"].value_counts()
world = world.merge(country_counts, left_on="name", right_index=True, how="left")
fig, ax = plt.subplots(1, 1, figsize=(20, 10))
world.plot(
    column="count",
    ax=ax,
    legend=True,
    missing_kwds={"color": "lightgrey"},
    cmap="coolwarm",
    legend_kwds={
        "label": "Number of genomes",
        "shrink": 0.5,
    },
    edgecolor="gray",
)
ax.set_title("Number of genomes per country")
ax.set_axis_off()
plt.tight_layout()
plt.savefig("figs/genomes_per_country.png", dpi=300)
plt.show()

# %%
