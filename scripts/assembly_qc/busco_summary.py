import argparse
import json
import pathlib
import re
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect summary statistics for different busco runs."
    )
    parser.add_argument("--fld", help="busco output folders", type=str, nargs="+")
    parser.add_argument("--csv", help="output csv file", type=str)
    return parser.parse_args()


def extract_info(fname):
    with open(fname, "r") as f:
        data = json.load(f)
    info = data["results"]
    return info


def extract_iso_name(fname):
    # Regular expression pattern
    pattern = r"enterobacterales_odb10\.(.*?)\.json$"

    # Search for the pattern
    match = re.search(pattern, fname)

    if match:
        return match.group(1)
    else:
        print("Warning: no match found for", fname)
        return None


if __name__ == "__main__":
    args = parse_args()

    # create summary dataframe
    df = []
    for fld in args.fld:
        # find json input file
        fld = pathlib.Path(fld)
        fnames = list(fld.glob("short_summary.specific.enterobacterales_odb10.*.json"))
        assert len(fnames) == 1
        fname = fnames[0]

        # extract isolate name
        iso = extract_iso_name(fname.name)
        # extract and append info
        info = extract_info(fname)
        info["iso"] = iso
        df.append(info)

    # create and save dataframe
    df = pd.DataFrame(df)
    df.to_csv(args.csv, index=False)
