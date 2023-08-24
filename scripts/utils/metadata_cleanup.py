# %%

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Clean up metadata")
    parser.add_argument("--metadata_in", type=str)
    parser.add_argument("--metadata_out", type=str)
    return parser.parse_args()


substitution_dict = {
    "collection_date": {"2015-01/2016-07": "2015"},
    "host": {
        "pantropical spotted dolphin": "dolphin",
        "Homo sapiens; 55 year old": "human",
        "Pig": "pig",
        "Homo sapiens": "human",
    },
    "geo_loc_name": {
        "USA: Minnesota": "USA",
        "USA: New York": "USA",
        "USA: Minneapolis": "USA",
        "USA: Burlington": "USA",
        "unknown": None,
        "Canada: Winnipeg": "Canada",
        "Brazil: Botucatu, Sao Paulo State": "Brazil",
        "USA:Utah": "USA",
        "USA:Houston": "USA",
        "United Kingdom:Manchester": "United Kingdom",
        "Australia:Brisbane": "Australia",
        "Czech Republic: Brno": "Czech Republic",
        "USA: San Diego": "USA",
        "USA: Houston": "USA",
        "China: Hainan Province": "China",
        "USA: Pittsburgh, Pennsylvania": "USA",
    },
    "isolation_source": {
        "wastewater treatment plant effluent": "wastewater",
        "missing": None,
        "unknown": None,
        "human, urinary tract infection": "UTI",
        "urine (asymptomatic UTI)": "UTI",
        "urine (recurrent UTI)": "UTI",
        "dog, faeces": "feces",
        "cecal sample": "feces",
        "vaginal swab": "UTI",
    },
}

if __name__ == "__main__":
    args = parse_args()

    df = pd.read_csv(args.metadata_in, index_col=0)
    for k, d in substitution_dict.items():
        df[k] = df[k].replace(d)
    df.to_csv(args.metadata_out)
# %%
