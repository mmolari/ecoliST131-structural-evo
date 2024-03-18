import pandas as pd
import numpy as np


def assign_category(df):
    df["cat"] = None
    keys = ["genomad", "integrons", "isescan"]

    for k, val in [
        ("IS", (False, False, True)),
        ("none", (False, False, False)),
        ("prophage", (True, False, True)),
        ("prophage", (True, False, False)),
        ("integron", (False, True, True)),
    ]:
        mask = np.all((df[keys] > 0) == val, axis=1)
        df.loc[mask, "cat"] = k
    return df
