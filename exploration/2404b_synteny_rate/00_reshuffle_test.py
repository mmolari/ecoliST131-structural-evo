# %%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
from Bio import Phylo

fname = "../../results/ST131_ABC/pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(fname, "newick")
ter_len = sum(l.branch_length for l in tree.get_terminals())
int_len = sum(l.branch_length for l in tree.get_nonterminals())
L_tot = ter_len + int_len
filt_aln_size = 2427416

print(f"internal length = {int_len}")
print(f"terminal length = {ter_len}")
print(f"total length = {L_tot}")
print(f"terminal / total ratio = {ter_len/L_tot}")

print(f"number of mutations = {filt_aln_size * L_tot}")
# %%
n_ter = 18
n_int = 3

# extra event
n_ter += 1
n_int += 1

# %%
# reshuffle test
N = n_ter + n_int
p_obs = n_ter / N
p_exp = ter_len / (ter_len + int_len)
bins = np.arange(N + 1)
b = sps.binom.cdf(bins, n=N, p=p_exp)

plt.plot(bins / N, b, "oC0")
plt.axvline(p_obs, color="k", ls="--", label="observed")
plt.axvline(p_exp, color="gray", ls="--", label="tree terminal fraction")
plt.legend()
plt.xlabel("fraction of terminal events")
plt.ylabel("cumulative probability")
plt.xlim(0.5, 1.05)
plt.tight_layout()
plt.show()

# %%
N = n_ter + n_int
L = ter_len + int_len
rate = N / L
ter_frac = ter_len / L

S = 100000
Nt = np.random.poisson(N * ter_frac, size=S)
Ni = np.random.poisson(N * (1 - ter_frac), size=S)
Ft = Nt / (Ni + Nt)

plt.hist(Ft, bins=1000, cumulative=True, density=True, histtype="step")
plt.axvline(p_obs, color="k", ls="--", label="observed")
plt.axvline(p_exp, color="gray", ls="--", label="tree terminal fraction")
plt.legend()
plt.xlabel("fraction of terminal events")
plt.ylabel("cumulative probability")
plt.xlim(0.5, 1.05)
plt.tight_layout()
plt.show()
# %%
