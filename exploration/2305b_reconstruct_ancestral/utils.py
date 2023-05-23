import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from collections import defaultdict
from itertools import combinations


def plot_venn(sets, keys, colors):
    fig, ax = plt.subplots(figsize=(6, 6))

    centers = dict(
        zip(keys, [np.array(xy) for xy in [(-1, 1), (1, 1), (1, -1), (-1, -1)]])
    )

    # plot four circles
    for i, k in enumerate(keys):
        ax.add_patch(
            mpl.patches.Circle(
                centers[k], 2, color=colors[k], alpha=0.3, zorder=0, label=keys[i]
            )
        )

    positions = {k: np.array(centers[k]) * 1.5 for k in keys}
    for k1, k2 in combinations(keys, 2):
        lab = "|".join(sorted([k1, k2]))
        positions[lab] = (centers[k1] + centers[k2]) / 1.4
    for k1, k2, k3 in combinations(keys, 3):
        lab = "|".join(sorted([k1, k2, k3]))
        positions[lab] = (centers[k1] + centers[k2] + centers[k3]) / 1.7
    positions["|".join(sorted(keys))] = np.array([0, 0])

    for i, (k, v) in enumerate(sets.items()):
        print(k, v)
        ax.text(*positions[k], v, ha="center", va="center", fontsize=12)

    for k in keys:
        p = centers[k] * 1.8
        ax.text(*p, k, ha="center", va="center", fontsize=12, color=colors[k])
    # plt.legend()
    plt.xlim(-3, 3)
    plt.ylim(-3, 3)
    for s in ax.spines:
        ax.spines[s].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    return fig
