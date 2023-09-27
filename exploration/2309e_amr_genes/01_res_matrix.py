# %%
import pandas as pd
import seaborn as sns
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
import os

# %%
dbase = "card"
# dbase = "ncbi"

# read the tree and extract isolate order
tree_file = "../../results/ST131/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
isolates = [leaf.name for leaf in tree.get_terminals()]

# Read the TSV file into a DataFrame
res_file = f"../../results/ST131/resistance/{dbase}_summary.txt"
df = pd.read_csv(res_file, sep="\t")

# Extract the relevant columns (excluding the first column 'FILE')
columns = df.columns[2:]  # Exclude the first two columns ('FILE' and 'NUM_FOUND')

# Replace '.' entries with 0
df[columns] = df[columns].replace(".", 0)


# Convert percentage strings to float numbers, and returns the number of them higher than a threshold


def parse_percentage_str(percentage_str, thr):
    percentages = percentage_str.split(";")
    percentages = [float(p) for p in percentages]
    return len([p for p in percentages if p >= thr])


threshold = 95
df[columns] = df[columns].applymap(
    lambda x: parse_percentage_str(x, threshold)
    if isinstance(x, str)
    else int(x > threshold)
)
col_order = df[columns].sum().sort_values(ascending=False).index

# Remove the '.gbk' suffix and path prefix from isolate names
df["#FILE"] = df["#FILE"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
df.set_index("#FILE", inplace=True)
df = df.loc[isolates]

# Set up the figure and axis for plotting
figsize = {
    "card": (15, 14),
    "ncbi": (10, 14),
}
grid_kw = {
    "card": {"width_ratios": [0.8, 3, 0.4]},
    "ncbi": {"width_ratios": [0.8, 3, 0.8]},
}
fig, axs = plt.subplots(
    1,
    3,
    figsize=figsize[dbase],
    gridspec_kw=grid_kw[dbase],
    sharey=False,
)

# Plot the tree without labels
ax = axs[0]
Phylo.draw(tree, axes=ax, show_confidence=False, label_func=lambda x: "", do_show=False)
ax.grid(alpha=0.5, axis="y")
ax.set_ylim(len(df) + 0.5, 0.5)
ax.set_ylabel("Isolates")

# Plot the matrix as a heatmap
ax = axs[1]
cax = ax.matshow(df[col_order], cmap="Blues")

# Customize the heatmap plot
ax.set_xticks(range(len(col_order)))
ax.set_yticks(range(len(df)))
ax.set_xticklabels(col_order, rotation=90)
ax.set_yticklabels(df.index)
ax.set_xlabel("resistance genes")
# ax.set_ylabel("Isolates")
plt.colorbar(cax, ax=axs[-1], label=f"n. copies (>{threshold}% identity)", shrink=0.5)
ax.set_ylim(len(df) - 0.5, -0.5)

# Plot the bar plot on the right side
ax = axs[2]
bar_positions = np.arange(len(df))
ax.barh(bar_positions, df["NUM_FOUND"], align="center", color="silver")
ax.set_yticks([])
ax.set_xlabel("total n. genes")
ax.set_ylim(len(df) - 0.5, -0.5)


# Adjust the layout to prevent overlap
sns.despine()
plt.tight_layout()

# make the right plot closer
plt.subplots_adjust(wspace=0.1)

# Display the plot
plt.savefig(f"figs/resistance_{dbase}.png", dpi=300)
# plt.savefig(f"figs/resistance_{dbase}.pdf")
plt.savefig(f"figs/resistance_{dbase}.svg")
plt.show()
# %%
