# %%
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

# very simple coalescent with selection for small genome size
beta = 0.1
N = 100
np.random.seed(0)
rand_new_genes = lambda: np.random.binomial(n=2, p=0.6) - 1


@dataclass
class Genome:
    n: int = 0  # number of extra genes / ISs

    def evolve(self):
        n_old = self.n
        self.n += rand_new_genes()
        self.n = max(0, self.n)
        delta = self.n - n_old
        return delta


def evolve(P):
    Ds = np.array([g.evolve() for g in P])
    N_mean = np.mean([g.n for g in P])
    E = np.exp(-beta * np.array([g.n - N_mean for g in P]))
    p = E / E.sum()
    P_idx = np.random.choice(range(N), size=N, p=p)
    # make a copy
    P_new = np.array([Genome(n=P[i].n) for i in P_idx])
    return P_new, Ds


P = np.array([Genome() for _ in range(N)])
T = 300
Gs = []
Ds = []
for t in range(T):
    P, D = evolve(P)
    G = np.array([g.n for g in P])
    Gs.append(G)
    Ds.append(D)

mean = np.array([G.mean() for G in Gs])
std = np.array([G.std() for G in Gs])

plt.plot(mean, label=r"$n(t)$")
plt.fill_between(range(T), mean - std, mean + std, alpha=0.5)

mean_d = np.array([D.mean() for D in Ds])
std_d = np.array([D.std() for D in Ds])

plt.plot(mean_d, label=r"$\Delta n(t)$")
plt.fill_between(range(T), mean_d - std_d, mean_d + std_d, alpha=0.5)
plt.axhline(0, color="black", linestyle="--")

plt.legend()
plt.xlabel("time")
plt.ylabel("number of ISs")
plt.tight_layout()
plt.show()

# %%
