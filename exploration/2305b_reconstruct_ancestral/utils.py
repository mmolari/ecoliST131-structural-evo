import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo

from collections import defaultdict
from itertools import combinations


def plot_tree_pairs(tree, pairs, iso_colors):
    # display tree with selected isolates in color
    fig, ax = plt.subplots(figsize=(6, 12))
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=lambda x: x.name if x.name in pairs else "",
        label_colors=lambda x: iso_colors[x] if x in pairs else "black",
        # branch_labels=lambda x: x.branch_length,
        # branch_labels_color="red",
        # branch_labels_size=10,
        show_confidence=False,
    )
    for s in ax.spines:
        ax.spines[s].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    return fig


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


def display_overlap(ctr, pairs, iso_colors):
    iso_position = {k: i for i, k in enumerate(pairs)}
    fig, axs = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

    # singleton blocks
    ax = axs[0]
    ct = [ctr[p] for p in pairs]
    bars = ax.bar(range(len(pairs)), ct, color=[iso_colors[i] for i in pairs])
    ax.bar_label(bars, ct)
    ax.set_title("only in one")

    # triplet blocks
    ax = axs[1]
    ct = []
    for p in pairs:
        triplet = [q for q in pairs if q != p]
        label = "|".join(sorted(triplet))
        ct.append(ctr[label])
    bars = ax.bar(range(len(pairs)), ct, color=[iso_colors[i] for i in pairs])
    ax.bar_label(bars, ct)
    ax.set_title("missing in one")

    # block pairs
    ax = axs[2]
    for k, v in ctr.items():
        isos = k.split("|")
        if len(isos) != 2:
            continue
        # draw arcs between isolate pairs
        iso_x = [iso_position[i] for i in isos]
        r = (np.max(iso_x) - np.min(iso_x)) / 2
        theta = np.arange(0, np.pi, np.pi / 100)
        xs = r * (np.cos(theta) - 1) + np.max(iso_x)
        ys = r * np.sin(theta)
        ax.plot(xs, ys, color="black")
        # add text on top of arc
        ax.text(np.mean(iso_x), r * 1.05, str(v), ha="center", va="bottom")
    ax.set_yticks([])
    ax.set_xticks(range(len(pairs)), pairs)
    ax.set_title("in pair")

    # in all
    all_label = "|".join(sorted(pairs))
    if all_label in ctr:
        n = ctr[all_label]
        ax.text(
            0.1,
            0.9,
            f"all = {n}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )

    # despine
    for ax in axs:
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)

    plt.tight_layout()
    return fig
