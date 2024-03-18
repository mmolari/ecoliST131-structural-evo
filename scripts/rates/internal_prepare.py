import argparse

import pandas as pd

import pypangraph as pp
import utils as ut


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--joints_dset", type=str)
    parser.add_argument("--joint_pangraphs", type=str, nargs="+")
    parser.add_argument("--out_AB_states", type=str)
    parser.add_argument("--out_event_info", type=str)
    return parser.parse_args()


def create_AB_states_df(jdf, pg_files):
    two_cat_isolates = {}

    AB_df = {}
    for j in jdf.index:

        pan = pp.Pangraph.load_json(pg_files[j])

        paths = ut.pangraph_to_path_dict(pan)
        path_cat = ut.path_categories(paths)

        assert len(path_cat) == 2, f"More than 2 categories in {j}"

        n1, p1, i1 = path_cat[0]
        n2, p2, i2 = path_cat[1]

        assert n2 > 1, "Only one path in category 2"

        two_cat_isolates[j] = (p1, p2)
        AB_df[j] = {i: "A" for i in i1} | {i: "B" for i in i2}
    AB_df = pd.DataFrame.from_dict(AB_df, orient="columns")
    AB_df.fillna(".", inplace=True)  # in case non-backbone edges are considered
    return AB_df, two_cat_isolates


def AB_nestedness(two_cat_isolates):
    """finds whether A<B, A>B (or possibly none of these A?B)"""
    nestedness = {}
    for j in two_cat_isolates:
        p1, p2 = two_cat_isolates[j]
        p1, p2 = set(p1), set(p2)
        if p1.issubset(p2):
            nestedness[j] = "A<B"
        elif p2.issubset(p1):
            nestedness[j] = "A>B"
        else:
            nestedness[j] = "A?B"
    nestedness = pd.Series(nestedness, name="event_type")
    return nestedness


if __name__ == "__main__":

    args = parse_args()

    # select non-terminal events
    jdf = pd.read_csv(args.joints_dset, index_col=0)
    mask = ~jdf["singleton"]
    jdf = jdf[mask]

    # pangraph joint file dictionary
    pg_files = ut.create_pangraph_file_dict(args.joint_pangraphs)

    # create PA states dataframe
    AB_df, two_cat_isolates = create_AB_states_df(jdf, pg_files)

    # infer gain/loss event type
    nestedness = AB_nestedness(two_cat_isolates)

    # save dataframes
    AB_df.to_csv(args.out_AB_states)
    nestedness.to_csv(args.out_event_info)

# %%
