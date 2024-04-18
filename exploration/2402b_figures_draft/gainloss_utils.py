import pandas as pd
import numpy as np
import pathlib


def assign_mge_category(df):
    df["cat"] = None
    keys = ["genomad", "integrons", "isescan", "defensefinder"]
    F, T = False, True
    for val, k in [
        ((T, F, T, F), "prophage"),
        ((T, F, F, F), "prophage"),
        ((T, F, F, T), "prophage"),
        ((F, T, T, T), "integron"),
        ((F, F, T, T), "defense"),
        ((F, F, T, F), "IS"),
        ((F, F, F, F), "none"),
    ]:
        mask = np.all((df[keys] > 0) == val, axis=1)
        df.loc[mask, "cat"] = k

    # check that no "NaN" is left
    assert df["cat"].isna().sum() == 0, "some categories are not assigned"

    # make ordered categorical variable
    df["cat"] = pd.Categorical(
        df["cat"],
        categories=["IS", "integron", "prophage", "defense", "none"],
        ordered=True,
    )

    return df


cat_colors = {
    "IS": "C0",
    "prophage": "C4",
    "integron": "C1",
    "none": "#b8b8b8",
    "defense": "C2",
}

# filenames
load_fld = pathlib.Path("../../results/ST131_ABC/rates/asm20-100-5")
fnames = {
    "cdf": {
        "terminal": load_fld / "terminal_coldspot.csv",
        "internal": load_fld / "nonsingleton_junct_df.csv",
    },
    "bdf": {
        "terminal": load_fld / "terminal_branches.csv",
        "internal": load_fld / "nonsingleton_branch_df.csv",
    },
    "coldspots": load_fld / "coldspot_df.csv",
    "events": load_fld / "merged_events_df.csv",
}


def load_cdf(fname):
    cdf = pd.read_csv(fname, index_col=0)
    cdf = assign_mge_category(cdf)
    return cdf
