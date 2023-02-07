import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="""merge distance dataframes""")
    parser.add_argument(
        "--dfs_in", type=str, nargs="+", help="input distance dataframes"
    )
    parser.add_argument("--csv", type=str, help="output distance dataframe")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    dfs = []
    for in_csv in args.dfs_in:
        df = pd.read_csv(in_csv).set_index(["si", "sj"])
        dfs.append(df)

    dfs = pd.concat(dfs, axis=1, verify_integrity=True, join="outer")
    dfs.to_csv(args.csv)
