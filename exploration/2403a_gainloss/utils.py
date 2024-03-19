import pandas as pd
import numpy as np


def assign_mge_category(df):
    df["cat"] = None
    keys = ["genomad", "integrons", "isescan", "defensefinder"]
    F, T = False, True
    for val, k in [
        ((F, F, T, F), "IS"),
        ((F, F, F, F), "none"),
        ((T, F, T, F), "prophage"),
        ((T, F, F, F), "prophage"),
        ((F, T, T, T), "integron"),
        ((F, F, T, T), "defense"),
        ((T, F, F, T), "prophage"),
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
