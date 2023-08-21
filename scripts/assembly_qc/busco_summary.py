import argparse
import json
import pathlib
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect summary statistics for different busco runs."
    )
    parser.add_argument("--fld", help="busco output folders", type=str, nargs="+")
    parser.add_argument("--csc", help="output csv file", type=str)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    df = []
    for fld in args.fld:
        fname = pathlib.Path(fld).glob("short_summary.specific.enterobacterales_odb10.*.json")

