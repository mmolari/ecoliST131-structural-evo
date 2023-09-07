# %%
import pandas as pd
import numpy as np
import scipy.stats as sps
import matplotlib.pyplot as plt
from functools import cache
from collections import defaultdict

df = pd.read_csv("data/event_df.csv")
T = df["branch_length"].to_numpy()
N = df["ev_count"].to_numpy()
# N = df["gain"].to_numpy()
# N = df["loss"].to_numpy()
plt.plot(T, N, ".", alpha=0.5)
plt.show()

lm = N.sum() / T.sum()


@cache
def log_factorial(n):
    if n == 0:
        return 0
    else:
        return log_factorial(n - 1) + np.log(n)


def logl(N, T, lm):
    lfN = [log_factorial(int(n)) for n in N]
    return sum(N * np.log(lm * T) - lm * T - lfN)


# %%
L_data = logl(N, T, lm)

LLs = []
for i in range(10000):
    N_sim = sps.poisson.rvs(lm * T)
    LLs.append(logl(N_sim, T, lm))
# %%

plt.hist(
    LLs, cumulative=True, bins=100, density=True, histtype="step", label="simulations"
)
plt.axvline(L_data, color="red", label="data")
plt.yscale("log")
plt.xlabel("log likelihood")
plt.ylabel("cumulative distribution")
plt.title("poisson log-likelyhood of number of gain events")
plt.tight_layout()
plt.legend()
plt.show()

# %%

LLs = {}
lm = {}
L_data = {}
for lab in ["ev_count", "gain", "loss"]:
    lm[lab] = N.sum() / T.sum()
    N = df[lab].to_numpy()
    L_data[lab] = logl(N, T, lm[lab])

    LLs[lab] = []
    for i in range(100000):
        N_sim = sps.poisson.rvs(lm[lab] * T)
        LLs[lab].append(logl(N_sim, T, lm[lab]))
# %%

fig, axs = plt.subplots(1, 3, figsize=(10, 3))
n = 0
for lab, txt in [("ev_count", "all"), ("gain", "gain"), ("loss", "loss")]:
    ax = axs[n]
    ax.hist(
        LLs[lab],
        cumulative=True,
        bins=100,
        density=True,
        histtype="step",
        label="simulations",
    )
    ax.axvline(L_data[lab], color="red", label="data")
    ax.set_yscale("log")
    ax.set_xlabel("poisson log-likelihood")
    ax.set_ylabel("cumulative distribution")
    ax.set_title(f"number of {txt} events")
    ax.text(
        0.1,
        0.9,
        r"$\lambda =$" + f"{lm[lab]:.3}",
        transform=ax.transAxes,
    )
    n += 1
ax.legend(loc="lower right")
plt.tight_layout()
plt.savefig(f"figs/log_likelihood.png")
plt.show()

# %%
