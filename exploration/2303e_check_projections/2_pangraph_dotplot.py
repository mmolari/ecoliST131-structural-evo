# %%
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns
import pypangraph as pp
from collections import defaultdict


def block_positions(pan, strain):
    p = pan.paths[strain]
    S = p.block_strands
    B = p.block_ids
    N = p.block_nums
    l = 0
    L = defaultdict(dict)
    for b, s, n in zip(B, S, N):
        lb = len(pan.blocks[b])
        L[b][n] = (l, l + lb) if s else (l + lb, l)
        l += lb
    return L


def block_dotplot(ax, pan, str1, str2, alpha_dupl=0.1):

    L1 = block_positions(pan, str1)
    L2 = block_positions(pan, str2)

    D1 = [b for b in L1.keys() if len(L1[b]) > 1]
    D2 = [b for b in L2.keys() if len(L2[b]) > 1]
    D = set(D1) | set(D2)

    p = pan.paths[str1]
    B1 = p.block_ids
    N1 = p.block_nums

    for b, n in zip(B1, N1):
        l1s, l1e = L1[b][n]
        c = "r" if b in D else "k"
        a = alpha_dupl if b in D else 1.0
        if b in L2:
            for (n, (l2s, l2e)) in L2[b].items():
                ax.plot([l1s, l1e], [l2s, l2e], color=c, alpha=a)
    ax.set_xlabel(str1)
    ax.set_ylabel(str2)


# %%

pangraph_file = "../../results/ST131/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(pangraph_file)

# %%

str_pairs = [
    ("NZ_JAOSCJ010000001", "NZ_CP076689"),
    ("NZ_CP104846", "NZ_CP104848"),
    ("NZ_CP035476", "NZ_CP035720"),
    ("NZ_CP103755", "NZ_CP049077"),
    ("NZ_SEVU01000007", "NZ_CP035476"),
    ("NZ_CP014497", "NZ_JAOSES010000001"),
    ("NZ_JAOSCG010000001", "NZ_CP021454"),
    ("NZ_CP014522", "NZ_JAOSCA010000001"),
    ("NZ_CP076687", "NZ_CP104846"),
    ("NZ_CP076687", "NZ_CP104848"),
    ("NZ_JAOSEJ010000001", "NZ_CP014497"),
]

for i1, i2 in str_pairs:

    if i1 < i2:
        i1, i2 = i2, i1

    fig_fld = pathlib.Path(f"figs/pairs/{i1}-{i2}-subplots")
    fig_fld.mkdir(exist_ok=True)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    block_dotplot(ax, pan, i1, i2, alpha_dupl=0.0)
    sns.despine(fig)
    plt.savefig(fig_fld / "nodupl_dotplot.png", dpi=300)
    plt.show()

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    block_dotplot(ax, pan, i1, i2, alpha_dupl=1.0)
    sns.despine(fig)
    plt.savefig(fig_fld / "dotplot.png", dpi=300)
    plt.show()
# %%
