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
        return "subset"
    elif b_in:
        return "start"
    elif e_in:
        return "end"
    else:
        return "no"


def within(B, E, b, e):
    if B <= E:
        return within_no_loop(B, E, b, e)
    else:
        return within_loop(B, E, b, e)
