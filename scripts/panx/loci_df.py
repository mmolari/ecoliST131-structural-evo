# %%
import argparse
import utils as ut

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_gbk", type=str, required=True)
    parser.add_argument("--out_df", type=str, required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    loci_df = ut.gbk_to_loci_df(args.in_gbk)
    loci_df.to_csv(args.out_df, index=False)
