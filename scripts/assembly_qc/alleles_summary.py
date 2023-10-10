import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="concatenate allele dataframes")
    parser.add_argument("--in_dfs", nargs="+", type=str)
    parser.add_argument("--summary_df", type=str)
    return parser.parse_args()


def load_allele(df_name):
    df = pd.read_csv(df_name, index_col=0, sep="\t")
    series = df["allele"].map(lambda x: str(int(x)) if pd.notnull(x) else x)
    series.name = df["locus"].iloc[0]
    return series


if __name__ == "__main__":
    args = parse_args()

    # load allele series
    srs = [load_allele(df_name) for df_name in args.in_dfs]

    # concatenate into dataframe
    df = pd.concat(srs, axis=1)

    df.to_csv(args.summary_df)
