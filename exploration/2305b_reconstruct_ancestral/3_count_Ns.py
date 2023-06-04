# %%
import os
from Bio import SeqIO
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# %%
accs = os.listdir("../../data/fa")
accs = [a.removesuffix(".fa") for a in accs]
# %%
ns = {}
for a in accs:
    seq = SeqIO.read(f"../../data/fa/{a}.fa", format="fasta").seq
    ns[a] = seq.count("N")
# %%
ns = pd.Series(ns)
# %%
ns[ns > 0]
