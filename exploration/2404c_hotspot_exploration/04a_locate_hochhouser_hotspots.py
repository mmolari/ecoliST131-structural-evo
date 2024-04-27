# %%
from Bio import SeqIO, Phylo
import pandas as pd
import matplotlib.pyplot as plt

fname = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(fname, "newick")
strains = [l.name for l in tree.get_terminals()]

fname = "../../results/ST131_ABC/pangraph/genome_lengths.csv"
df = pd.read_csv(fname, index_col=0)

genes = ["nanS", "yjgB", "yjgR"]
colors = ["red", "green", "blue"]


def find_gene(genes, fname):
    features = []
    for rec in SeqIO.parse(fname, "genbank"):
        for feature in rec.features:
            if feature.type == "CDS":
                for gene in genes:
                    if gene in feature.qualifiers.get("gene", []):
                        features.append(feature)
    return features


locations = {}
for iso, row in df.iterrows():
    fname = f"../../data/gbk/{iso}.gbk"
    feat = find_gene(genes, fname)
    locations[iso] = feat

# %%
gene_col = dict(zip(genes, colors))

fig, ax = plt.subplots(1, 1, figsize=(8, 12))
y = 0
for iso in strains:
    L = df.loc[iso, "length"]
    y -= 1
    center = None
    for feat in locations[iso]:
        if feat.qualifiers["gene"][0] == genes[0]:
            strand = feat.location.strand
            center = feat.location.start if (strand == 1) else feat.location.end
    if center is None:
        continue
    for feat in locations[iso]:
        start = feat.location.start
        end = feat.location.end
        if strand == 1:
            start -= center
            end -= center
        else:
            start = center - start
            end = center - end
        start %= L
        end %= L
        ax.plot([start, end], [y, y], "+", color=gene_col[feat.qualifiers["gene"][0]])
plt.tight_layout()
plt.show()

# %%

# I think I found the spot:
# RYYAQMEJGY_r__ZTHKZYHPIX_f
# hochhouser record: CP014348.1
