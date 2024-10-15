# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import pathlib
from collections import defaultdict
import pypangraph as pp

# look into inversion that happened multiple times
# 115 kbps
# isolates:
# NZ_CP107114.1
# NZ_CP107184.1
# NZ_CP107182.1

fld = pathlib.Path("../../results/ST131_ABC")
fig_fld = pathlib.Path("figs/f00")
fig_fld.mkdir(exist_ok=True, parents=True)

colors_file = fld / "pangraph/coresynt-asm20-100-5/blocks.csv"
mergers_file = fld / "pangraph/coresynt-asm20-100-5/mergers.csv"

cl_df = pd.read_csv(colors_file, index_col=0)
mg_df = pd.read_csv(mergers_file, index_col=0)

pangraph_file = fld / "pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()
# %%
iso_A = "NZ_CP107114.1"
iso_B = "NZ_CP107184.1"
iso_C = "NZ_CP107182.1"


def position_dict(path):
    pos = path.block_positions
    blocks = path.block_ids
    strand = path.block_strands
    pos_dict = defaultdict(list)
    for i, b in enumerate(blocks):
        pos_dict[b].append((pos[i], pos[i + 1], strand[i]))
    return pos_dict


def dotplot(pan, isoA, isoB, bdf, cl_df, mg_df):
    pA = pan.paths[isoA]
    pB = pan.paths[isoB]

    posA = position_dict(pA)
    posB = position_dict(pB)

    fig, ax = plt.subplots(figsize=(10, 10))
    for b, PA in posA.items():
        if b not in posB:
            continue
        for bA, eA, sA in PA:
            for bB, eB, sB in posB[b]:
                if (bA - eA) > 2e6 or (bB - eB) > 2e6:
                    continue
                x = [bA, eA]
                y = [bB, eB] if sA == sB else [eB, bB]
                col = "gray"
                lw = 0.5
                if bdf.loc[b, "core"] and bdf.loc[b, "len"] > 500:
                    if b in mg_df.index:
                        mg = mg_df.loc[b, "0"]
                    else:
                        mg = b
                    col = cl_df.loc[mg, "color"]
                    lw = 1
                    if mg == "RYYAQMEJGY":
                        lw = 3
                ax.plot(x, y, color=col, lw=lw)
    ax.set_xlabel(isoA)
    ax.set_ylabel(isoB)
    ax.set_aspect("equal")
    plt.grid(alpha=0.3)

    return fig, ax


fig, ax = dotplot(pan, iso_A, iso_B, bdf, cl_df, mg_df)
plt.tight_layout()
plt.savefig(fig_fld / "dotplot_114_184_full.png")
plt.show()

# %%

fig, ax = dotplot(pan, iso_A, iso_B, bdf, cl_df, mg_df)
i = 5e4
j = -1e5
ax.set_xlim(2.4e6 + i, 3e6 + i)
ax.set_ylim(4.4e6 + j, 5e6 + j)
plt.tight_layout()
plt.savefig(fig_fld / "dotplot_114_184_zoom.png")
plt.show()

# %%

fig, ax = dotplot(pan, iso_A, iso_C, bdf, cl_df, mg_df)
plt.tight_layout()
plt.savefig(fig_fld / "dotplot_114_182_full.png")
plt.show()

# %%

fig, ax = dotplot(pan, iso_A, iso_C, bdf, cl_df, mg_df)
i = 0
j = 0
ax.set_xlim(2.4e6 + i, 3e6 + i)
ax.set_ylim(0.2e6 + j, 0.8e6 + j)
plt.tight_layout()
plt.savefig(fig_fld / "dotplot_114_182_zoom.png")
plt.show()

# %%

fig, ax = dotplot(pan, iso_B, iso_C, bdf, cl_df, mg_df)
plt.tight_layout()
plt.savefig(fig_fld / "dotplot_184_182_full.png")
plt.show()


# %%

fig, ax = dotplot(pan, iso_B, iso_C, bdf, cl_df, mg_df)
i = 0
j = 0
ax.set_xlim(4.2e6 + i, 5e6 + i)
ax.set_ylim(0e6 + j, 0.8e6 + j)
plt.tight_layout()
plt.savefig(fig_fld / "dotplot_184_182_zoom.png")
plt.show()

# %%

fig, ax = dotplot(pan, iso_A, iso_B, bdf, cl_df, mg_df)
ax.set_xlim(2.62e6, 2.64e6)
ax.set_ylim(4.42e6, 4.44e6)
plt.tight_layout()
plt.show()

# %%
