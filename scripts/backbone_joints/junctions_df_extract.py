import pandas as pd
import pypangraph as pp
import argparse
import utils as ut


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Given a pangraph, generates dataframes with information on core-edge
        length and frequency."""
    )
    parser.add_argument("--pangraph", type=str, required=True)
    parser.add_argument("--len_thr", type=int, required=True)
    parser.add_argument("--df_len", type=str, required=True)
    parser.add_argument("--df_freq", type=str, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    pan = pp.Pangraph.load_json(args.pangraph)
    bdf = pan.to_blockstats_df()
    paths = ut.pangraph_to_path_dict(pan)

    def is_core(node_id):
        return (bdf.loc[node_id, "len"] >= args.len_thr) and bdf.loc[node_id, "core"]

    # core-junctions dataframe
    jdf = []
    for iso, path in paths.items():
        junctions = ut.path_junction_split(path, is_core)
        for J in junctions:
            edge = J.flanking_edge()
            L = sum(bdf.loc[node.id, "len"] for node in J.center.nodes)
            jdf.append({"iso": iso, "edge": edge.to_str_id(), "len": L})
    jdf = pd.DataFrame(jdf)
    jdf = jdf.pivot_table(index="iso", columns="edge", values="len")

    # extract edge frequency
    edge_freq = jdf.notna().sum(axis=0).sort_values(ascending=False)
    edge_freq.name = "count"
    edge_freq.to_csv(args.df_freq)

    # sort by edge frequency and save
    jdf = jdf[edge_freq.index]
    jdf.to_csv(args.df_len)
