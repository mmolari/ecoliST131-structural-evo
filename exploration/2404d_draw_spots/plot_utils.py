import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio import SeqIO, Phylo
import pypangraph as pp
from collections import defaultdict


def load_tree():
    fname = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
    tree = Phylo.read(fname, "newick")
    strains = [l.name for l in tree.get_terminals()]
    return tree, strains


def load_coldspots_df():
    fname = "../../results/ST131_ABC/rates/asm20-100-5/coldspot_df.csv"
    df = pd.read_csv(fname)
    return df


def load_pangraph(edge):
    fname = "../../results/ST131_ABC/backbone_joints/asm20-100-5/joints_pangraph/"
    fname += f"{edge}.json"
    pan = pp.Pangraph.load_json(fname)
    return pan


def load_MGEs(edge):
    dfs = []
    for mge, kind in [
        ("ISEScan", "IS"),
        ("genomad", "prophage"),
        ("defensefinder", "defense system"),
        ("integronfinder", "integron"),
    ]:
        fname = (
            f"../../results/ST131_ABC/annotations/junct_pos_asm20-100-5/{mge}_real.csv"
        )
        df = pd.read_csv(fname)
        mask = df["junction"] == edge
        df["MGE_type"] = kind
        dfs.append(df[mask])
    return pd.concat(dfs, axis=0)
