# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pypangraph as pp
import utils as ut

len_thr = 500
dset = "ST131_full"
edge_file = f"../../results/{dset}/backbone_joints/asm20-100-5/core_edges.csv"
pangraph_file = f"../../results/{dset}/pangraph/asm20-100-5-polished.json"

# %%
pan = pp.Pangraph.load_json(pangraph_file)
bdf = pan.to_blockstats_df()
edges_df = pd.read_csv(edge_file)

# %%
paths = ut.pangraph_to_path_dict(pan)

# %%


def encompassing_edges(path, core_f):
    res = {}

    first_bid = None
    done = False
    L = len(path)
    for i, node in enumerate(path.nodes):
        if core_f(node.id):
            first_bid = node.id
            break
    edge = ut.Edge(path.nodes[i], None)
    link = []
    i += 1
    while not done:
        node = path.nodes[i]
        is_core = core_f(node.id)
        if is_core:
            edge.right = node
            res[edge] = link
            link = []
            edge = ut.Edge(node, None)
            if node.id == first_bid:
                done = True
        else:
            link.append(node)
        i += 1
        i %= L
    return res


def core_f(block_id):
    keep = bdf.loc[block_id, "core"]
    keep &= bdf.loc[block_id, "len"] > len_thr
    return keep


links = {}
for iso, path in paths.items():
    links_bid = encompassing_edges(path, core_f)
    links_len = {
        e: sum([bdf.loc[node.id, "len"] for node in l]) for e, l in links_bid.items()
    }
    links[iso] = links_len
# %%

df = []
for iso, link in links.items():
    for edge, L in link.items():
        res = {}
        res["edge"] = edge
        res["iso"] = iso
        res["len"] = L
        df.append(res)
df = pd.DataFrame(df)
df = df.pivot(index="iso", columns="edge", values="len")
df.columns = [e.to_str_id() for e in df.columns]
col_ord = df.notna().sum().sort_values(ascending=False).index
df = df[col_ord]
df

# %%
fig, ax = plt.subplots(figsize=(20, 20))
sns.heatmap(df.notna(), cmap="Blues", ax=ax)
plt.show()


fig, ax = plt.subplots(figsize=(20, 20))
sns.heatmap(df, cmap="rainbow", ax=ax)
plt.show()
# %%
core_edges_1 = edges_df[edges_df["count"] == len(pan.strains())]["edge"].to_numpy()
core_edges_2 = df.columns[df.notna().sum() == len(pan.strains())].to_numpy()

ce1 = set([ut.Edge.from_str_id(e) for e in core_edges_1])
ce2 = set([ut.Edge.from_str_id(e) for e in core_edges_2])

assert ce1 == ce2
# %%
