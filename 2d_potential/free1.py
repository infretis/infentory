import MDAnalysis as mda
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

fig, axs = plt.subplots(2, 3, figsize=(10.5, 2.625*2), layout="constrained")

for idx, ax in enumerate(axs[0, :]):
    ax.set_xlim([-0.5, 0.5])
    ax.set_xticklabels([])
    if idx > 0:
        ax.set_yticklabels([])
for idx, ax in enumerate(axs[1, :]):
    ax.set_xlim([-0.5, 0.5])
    # ax.set_xticklabels([])
    if idx > 0:
        ax.set_yticklabels([])
axs[0, 0].pcolormesh(x, y, -np.log(h1), rasterized=True)
axs[0, 1].pcolormesh(x, y, -np.log(h3), rasterized=True)
axs[0, 2].pcolormesh(x, y, -np.log(h2), rasterized=True)

if os.path.exists("wham2d/histo_committor.txt"):
    axs[1, 0].pcolormesh(x, y, c1, rasterized=True)
    axs[1, 1].pcolormesh(x, y, c3, rasterized=True)
    axs[1, 2].pcolormesh(x, y, c2, rasterized=True)


axs[1, 0].set_title("Conditional Committor Surface")
axs[1, 1].set_title("Committor Surface")
axs[1, 2].set_title("Conditional Committor Surface")
axs[0, 0].set_title("Conditional Free Energy")
axs[0, 1].set_title("Free Energy")
axs[0, 2].set_title("Conditional Free Energy")

axs[1, 0].set_xlabel("x")
axs[1, 1].set_xlabel("x")
axs[1, 2].set_xlabel("x")

axs[0, 0].set_ylabel("y")
axs[1, 0].set_ylabel("y")

plt.savefig("2d1.pdf", dpi=600)
plt.show()
