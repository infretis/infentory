import matplotlib.pyplot as plt
import numpy as np
import os
x = np.loadtxt("wham2d/histo_xval.txt")
y = np.loadtxt("wham2d/histo_yval.txt")
h1 = np.loadtxt("wham2d/histo_probability.txt").T
h2 = h1.copy()
h2 = h2[::-1,::-1]

h3 = np.array([h1, h2])
masked_data = np.ma.masked_array(h3, np.isinf(h3))
h3 = np.ma.average(masked_data, axis=0)

if os.path.exists("wham2d/histo_committor.txt"):
    c1 = np.loadtxt("wham2d/histo_committor.txt").T
    c2 = c1.copy()
    c2 -= 1
    c2[c2!=np.inf] = np.abs(c2[c2!=np.inf])
    c2 = c2[::-1,::-1]

    c3 = np.array([c1, c2])
    masked_data = np.ma.masked_array(c3, np.isinf(c3))
    c3 = np.ma.average(masked_data, axis=0)

# fig, axs = plt.subplots(4, 3, figsize=(10.5, 2.625*4), layout="constrained")
fig, axs = plt.subplots(4, 3, layout="constrained")

for idx, ax in enumerate(axs[0, :]):
    ax.set_xlim([-1.3, 1.3])
    ax.set_ylim([-1.3, 20])
    ax.set_xticklabels([])
    if idx > 0:
        ax.set_yticklabels([])

for idx, ax in enumerate(axs[1, :]):
    ax.set_xlim([-1.3, 1.3])
    ax.set_xticklabels([])
    if idx > 0:
        ax.set_yticklabels([])
for idx, ax in enumerate(axs[2, :]):
    ax.set_xlim([-1.3, 1.3])
    ax.set_xticklabels([])
    if idx > 0:
        ax.set_yticklabels([])
for idx, ax in enumerate(axs[3, :]):
    ax.set_xlim([-1.3, 1.3])
    ax.axvline(0, color="k")
    ax.axhline(0.5, color="k")
    if idx > 0:
        ax.set_yticklabels([])

# x1d = np.loadtxt("../wham1d/histo_xval.txt")
xd = np.loadtxt("wham1d/histo_xval.txt")
free1d = np.loadtxt("wham1d/histo_free_energy.txt")
free1d2 = np.loadtxt("wham1d/histo_probability.txt")

if os.path.exists("wham2d/histo_committor.txt"):
    commit1 = np.loadtxt("wham1d/histo_committor.txt")
free1d3 = (free1d2 + free1d2[::-1])/2
free1d4 = -np.log(free1d3)

print(xd)

def pot(a, b, c, x):
    return a*x**4 - b*(x-c)**2

y1 = pot(a=1, b=2, c=0, x=xd)
axs[0, 0].plot(xd, (y1-np.min(y1))*14.285714285714285, color="k", alpha=0.8, ls="--", lw=2, label="analytical")
axs[0, 1].plot(xd, (y1-np.min(y1))*14.285714285714285, color="k", alpha=0.8, ls="--", lw=2)
axs[0, 2].plot(xd, (y1-np.min(y1))*14.285714285714285, color="k", alpha=0.8, ls="--", lw=2)

axs[0, 0].legend(loc="upper left")
axs[0, 0].plot(xd, free1d, lw=2)
axs[0, 1].plot(xd, free1d4 - np.min(free1d4), lw=2)
axs[0, 2].plot(xd[::-1], free1d, lw=2)

axs[1, 0].pcolormesh(x, y, -np.log(h1), rasterized=True)
axs[1, 1].pcolormesh(x, y, -np.log(h3), rasterized=True)
axs[1, 2].pcolormesh(x, y, -np.log(h2), rasterized=True)

if os.path.exists("wham2d/histo_committor.txt"):
    axs[2, 0].pcolormesh(x, y, c1, rasterized=True)
    axs[2, 1].pcolormesh(x, y, c3, rasterized=True)
    axs[2, 2].pcolormesh(x, y, c2, rasterized=True)
    axs[3, 0].plot(xd, commit1)
    axs[3, 2].plot(xd, np.abs(1-commit1[::-1]))
    axs[3, 1].plot(xd, (commit1+np.abs(1-commit1[::-1]))/2)


axs[0, 0].set_title("Conditional")
axs[0, 1].set_title("Unconditional")
axs[0, 2].set_title("Conditional")
axs[2, 0].set_title("Conditional Committor Surface")
axs[2, 1].set_title("Committor Surface")
axs[2, 2].set_title("Conditional Committor Surface")
axs[1, 0].set_title("Conditional Free Energy")
axs[1, 1].set_title("Free Energy")
axs[1, 2].set_title("Conditional Free Energy")

axs[0, 0].set_ylabel("(Conditional) Free Energy")
axs[1, 0].set_ylabel("Velocity")
axs[2, 0].set_ylabel("Velocity")
axs[3, 0].set_ylabel("(Conditional) Committor")

axs[3, 0].set_xlabel("x")
axs[3, 1].set_xlabel("x")
axs[3, 2].set_xlabel("x")

# plt.savefig("double_well.pdf", dpi=600)
plt.show()
