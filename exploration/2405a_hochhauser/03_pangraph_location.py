# %%
import pypangraph as pp
import pandas as pd
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib

fig_fld = pathlib.Path("figs/n03")
fig_fld.mkdir(exist_ok=True, parents=True)

fname = "../../results/ST131_ABC/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(fname)
bdf = pan.to_blockstats_df()
loc = pp.Locator(pan)

fname = "res/flanking_genes_loc.csv"
df = pd.read_csv(fname)

strain = "NZ_CP096110.1"
# %%
hs_locs = []
locations = defaultdict(dict)
for i, gr in df.groupby("hotspot"):
    if len(gr) != 2:
        continue
    g1, g2 = gr.iloc[0], gr.iloc[1]
    s1, e1 = g1["start"], g1["end"]
    s2, e2 = g2["start"], g2["end"]
    E = max(e1, e2)
    S = min(s1, s2)
    res1 = loc.find_interval(strain, s1, e1)
    res2 = loc.find_interval(strain, s2, e2)
    res3 = loc.find_interval(strain, S, E)
    locations[i][g1.coregene] = res1
    locations[i][g2.coregene] = res2
    locations[i]["all"] = res3
    assert np.abs(E - S) < 1e6
    hs_locs.append({"id": i, "iso": strain, "beg": S, "end": E, "type": "hotspot"})

hs_locs = pd.DataFrame(hs_locs)
hs_locs.to_csv("res/hs_genome_loc.csv", index=False)
# %%
res = []
coreintervals = defaultdict(list)
for i, locs in locations.items():
    for k, v in locs.items():
        if k == "all":
            for b, c, o in zip(v[0], v[1], v[2]):
                is_core = bdf.loc[b, "core"]
                is_long = bdf.loc[b, "len"] >= 500
                if is_core and is_long:
                    coreintervals[i].append((b, o[2]))

        else:
            B, C, O = v
            for b, c, o in zip(B, C, O):
                res.append(
                    {
                        "gene": k,
                        "hotspot": i,
                        "block": b,
                        "core": bdf.loc[b, "core"],
                        "n. strains": bdf.loc[b, "n. strains"],
                    }
                )
res = pd.DataFrame(res)
res
# %%
fig, ax = plt.subplots(1, 1, figsize=(3, 7))

for (hs, ud), gr in res.groupby(["hotspot", "gene"]):
    x = 0
    y = hs
    x = 0 if ud == "upstream" else 1
    core = np.any(gr["core"])
    nstrains = np.min(gr["n. strains"])
    marker = "s"
    color = "green" if core else "black"
    ax.scatter(x, y, color=color, marker=marker, s=80)
    if not core:
        ax.annotate(
            f"{nstrains:.0f}",
            (x, y),
            textcoords="offset points",
            xytext=(10, 0),
            va="center",
        )

ax.scatter([], [], color="green", marker="s", label="in core block(s)")
ax.scatter([], [], color="black", marker="s", label="in non-core block(s)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.05))


ax.set_xticks([0, 1])
ax.set_xticklabels(["upstream", "downstream"])
ax.set_xlim(-0.5, 1.5)
ax.set_title("Flanking genes")
sns.despine()
plt.tight_layout()
plt.savefig(fig_fld / "flanking_genes_core.png", dpi=200)
plt.show()

# %%
