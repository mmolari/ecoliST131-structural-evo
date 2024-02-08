import json
import pandas as pd
import numpy as np
import argparse

#       B                        E
#       |------------------------|
#            b           e
#            |-----------|                   subset
#     |-----------------------------|        superset
#                       |----------------|   start
#    |--------------|                        end
#                                     |----| no


def within_no_loop(B, E, b, e):
    """Given an target interval [B, E], with B <= E and a query interval [b, e], with b < e,
    determines whether the query interval overlaps with the target."""
    assert B <= E, f"{B=} > {E=}"
    assert b <= e, f"{b=} > {e=}"
    if b >= B and e <= E:
        return "subset"
    elif b <= B and e >= E:
        return "superset"
    elif b >= B and b <= E:
        return "start"
    elif e >= B and e <= E:
        return "end"
    else:
        return "no"


#   0   E                        B    L
#   ----|                        |-----
#            b           e
#            |-----------|                   no
#   |-----|                                  start
#                               |-----|      end
#   |-|                                      subset


def within_loop(B, E, b, e):
    """Given an target interval [B, L] union [0, E], with B > E and a query interval [b, e],
    with b < e, determines whether the query interval overlaps with the target."""
    assert B > E, f"{B=} < {E=}"
    assert b <= e, f"{b=} > {e=}"
    b_in = (b >= B) or (b <= E)
    e_in = (e >= B) or (e <= E)
    if b_in and e_in:
        if b < E and e > B:
            raise AssertionError(f"within loop: wrap-around {b=} {e=} {B=} {E=}")
        return "subset"
    elif b_in:
        return "start"
    elif e_in:
        return "end"
    else:
        return "no"


def rotate(x, delta, L):
    return (x + delta) % L


def within(B, E, b, e, L):

    if b > e:
        delta = -b
        B, E, b, e = (
            rotate(B, delta, L),
            rotate(E, delta, L),
            rotate(b, delta, L),
            rotate(e, delta, L),
        )

    if B <= E:
        return within_no_loop(B, E, b, e)
    else:
        return within_loop(B, E, b, e)


def tests():
    L = 100
    wi = lambda B, E, b, e: within(B, E, b, e, L)

    #       B              E            L
    #       |==============|
    #              |--|                    subset
    #              |------------|          start
    # --|          |---------------------  start
    #    |-------------|                   end
    # --------|                   |------  end
    #   |-----------------------------|    superset
    # ---------------------------|    |--  superset
    # --| |------------------------------  superset
    #                         |--------|   no
    # --|                           |----  no
    # --------|          |---------------  ??

    assert wi(30, 60, 35, 45) == "subset"
    assert wi(30, 60, 35, 65) == "start"
    assert wi(30, 60, 35, 10) == "start"
    assert wi(30, 60, 10, 45) == "end"
    assert wi(30, 60, 90, 45) == "end"
    assert wi(30, 60, 25, 65) == "superset"
    assert wi(30, 60, 90, 70) == "superset"
    assert wi(30, 60, 25, 10) == "superset"
    assert wi(30, 60, 80, 85) == "no"
    assert wi(30, 60, 80, 10) == "no"
    try:
        wi(30, 60, 55, 35)
    except AssertionError:
        pass
    else:
        assert False

    # 0   E                        B    L
    # ====|                        |=====
    #                                |--|  subset
    # -|                             |---  subset
    # |--|                                 subset
    # |-------|                            start
    # ------------------|            |---  start
    #                           |-----|    end
    # --|                       |--------  end
    #        |-------------|               no
    #  |-----------------------------|     ??

    assert wi(60, 30, 65, 90) == "subset"
    assert wi(60, 30, 65, 10) == "subset"
    assert wi(60, 30, 10, 20) == "subset"
    assert wi(60, 30, 10, 40) == "start"
    assert wi(60, 30, 90, 40) == "start"
    assert wi(60, 30, 35, 70) == "end"
    assert wi(60, 30, 35, 10) == "end"
    assert wi(60, 30, 40, 50) == "no"
    try:
        wi(60, 30, 10, 90)
    except AssertionError:
        pass
    else:
        assert False


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--iso_len", help="Input file with isolate lengths")
    parser.add_argument(
        "--junction_pos_json", help="Input file with junction positions"
    )
    parser.add_argument("--element_pos_df", help="Input file with element positions")
    parser.add_argument(
        "--output_pos", help="Output file with element junction assignment"
    )
    parser.add_argument("--random", action="store_true", help="Use random positions")
    parser.add_argument(
        "--zero-based-input-pos",
        action="store_true",
        help="Input positions are 0-based",
    )
    return parser.parse_args()


# assumptions:
# - junction positions are 0-based
# - element positions are 1-based unless zero-based-input-pos is set

if __name__ == "__main__":

    tests()

    args = parse_args()

    # load IS dataframe
    df_element = pd.read_csv(args.element_pos_df, index_col=0)

    # load joint coordinates dictionary
    with open(args.junction_pos_json) as f:
        junction_positions = json.load(f)

    # isolates genome length
    iso_L = pd.read_csv(args.iso_len, index_col=0)["length"].to_dict()

    rand_pos = args.random
    np.random.seed(42)

    zero_based = args.zero_based_input_pos

    results = []

    for el_id, row in df_element.iterrows():
        B, E = row["beg"], row["end"]
        iso = row["iso"]
        L = iso_L[iso]

        # make 0-based
        if not zero_based:
            B, E = B - 1, E - 1

        # optional randomize position
        if rand_pos:
            delta_L = np.random.randint(L)
            B, E = (B + delta_L) % L, (E + delta_L) % L

        # locate on joints
        for j, pos_d in junction_positions.items():

            if not iso in pos_d:
                continue

            cb, ab, ae, ce, strand = pos_d[iso]

            in_core = within(cb, ce, B, E, L)
            if in_core == "no":
                continue
            in_acc = within(ab, ae, B, E, L)
            if in_acc == "no":
                continue
            results.append(
                {
                    "id": el_id,
                    "junction": j,
                    "iso": iso,
                    "in_core": in_core,
                    "in_acc": in_acc,
                    "ib": B,
                    "ie": E,
                    "jcb": cb,
                    "jce": ce,
                    "jab": ab,
                    "jae": ae,
                    "j_strand": strand,
                }
            )

    results = pd.DataFrame(results)
    results.to_csv(args.output_pos, index=False)
