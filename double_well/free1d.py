import matplotlib.pyplot as plt
import numpy as np
#plt.savefig('dump.pdf')
#from inftools.misc.data_helper import data_reader
#from inftools.misc.infinit_helper import read_toml
import os


fig, axs = plt.subplots(1, 3, figsize=(10.5, 2.625*1), layout="constrained")


xd = np.loadtxt("wham1d/histo_xval.txt")
free1d = np.loadtxt("wham1d/histo_free_energy.txt")
free1d2 = np.loadtxt("wham1d/histo_probability.txt")
if os.path.exists("wham1d/histo_committor.txt"):
    commit1 = np.loadtxt("wham1d/histo_committor.txt")

pcross = np.loadtxt("wham1d/Pcross.txt")
free1d3 = (free1d2 + free1d2[::-1])/2
free1d4 = -np.log(free1d3)

def pot(a, b, c, x):
    return a*x**4 - b*(x-c)**2

y1 = pot(a=1, b=2, c=0, x=xd)

axs[0].set_xlim([-1.3, 1.3])
axs[0].set_ylim([-1.3, 20])

axs[0].set_xlabel("x")
axs[1].set_xlabel("x")
axs[2].set_xlabel("x")
axs[0].set_ylabel("(Conditional) Free Energy [kbT]")
axs[1].set_ylabel("(Conditional) Committor")
axs[2].set_ylabel("Conditional Crossing Probabillity")


axs[1].set_xlim([-1.3, 1.3])
axs[1].axvline(0, color="k", alpha=0.2)
axs[1].axhline(0.5, color="k", alpha=0.2)

axs[2].plot(pcross[:, 0],    pcross[:, -1], color="C0")
axs[2].plot(pcross[::-1, 0], pcross[:, -1], color="C1")
axs[2].set_yscale("log")

def dw(a, b, c, x):
    return a*x**4 - b*(x-c)**2

y1 = dw(a=1, b=2, c=0, x=xd)

# axs[0].legend(loc="upper left")
axs[0].plot(xd, free1d, lw=2, label=r"cond $\rightarrow$", color="C0", alpha=0.5)
axs[0].plot(xd[::-1], free1d, lw=2, label=r"cond $\leftarrow$", color="C1", alpha=0.5)
axs[0].plot(xd, free1d4 - np.min(free1d4), lw=2, label="Uncond", color="r")
axs[0].plot(xd, (y1-np.min(y1))*14.285714285714285, color="k", alpha=0.8, ls="--", lw=2, label="analytical")

if os.path.exists("wham1d/histo_committor.txt"):
    axs[1].plot(xd, commit1, color="C0", alpha=0.5)
    axs[1].plot(xd, np.abs(1-commit1[::-1]), color="C1", alpha=0.5)
    axs[1].plot(xd, (commit1+np.abs(1-commit1[::-1]))/2, color="r")
# axs[0].legend()
fig.legend(
    *axs[0].get_legend_handles_labels(),
    loc="upper center",
    bbox_to_anchor=(0.5, 1.11),
    ncol=4,           # number of legend columns
    frameon=False
)

# plt.savefig("1d.pdf")
plt.show()
