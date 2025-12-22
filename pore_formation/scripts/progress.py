import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import os

from inftools.misc.infinit_helper import read_toml


fig, axes = plt.subplots(1, 2, figsize=(7.0, (1.800)*2), constrained_layout=True)

axes[1].set_yticklabels([])
axes[0].set_ylim([1, 3.5])
axes[0].set_xlim([0.2, 2.2])
axes[1].set_ylim([1, 3.5])
axes[1].set_xlim([0.2, 2.2])
axes[0].set_ylabel(r"$d_{\rm LLP}$ [nm]")
axes[0].set_xlabel(r"$\xi_{\rm p}$")
axes[1].set_xlabel(r"$\xi_{\rm p}$")


x = np.loadtxt("y.txt")
y = np.loadtxt("x.txt")/10
ufree = np.loadtxt("ufree.txt")
ufree = ufree - np.min(ufree)
commi = np.loadtxt("committor.txt")

a = 0.7577519583264714
x0 = np.linspace(0, 3, 100)
c = [-2.3, 0.8]
for c0 in c:
    y0 = c0 - (1/a)*x0
    axes[0].plot(x0, -y0, ls="--", zorder=3000, color="k")
    axes[1].plot(x0, -y0, ls="--", zorder=3000, color="k")

axes[0].contour(x, y, ufree, colors='black', alpha=0.25, zorder=2000, levels=20)
axes[1].contour(x, y, commi, colors='black', alpha=0.25, zorder=2000, levels=20)
fifty_mask = (0.40 < commi) & (commi < 0.60)
fifty_mask &= np.isfinite(commi)
avg_hist_fifty = commi.copy()
avg_hist_fifty[~fifty_mask] = np.inf
axes[1].contourf(y, x, avg_hist_fifty, colors='r', alpha=0.25, zorder=100, levels=20)


x, y, xs, ys = [], [], [], []
reactive = False
for i in range(0, 1000):
    path = f"../load/{i}/order.txt"
    if not os.path.exists(path):
        break
    data = np.loadtxt(path)
    if data[-1, 1] > 0:
        reactive = True
        break
    x += list(data[:, 2]) + [None]
    y += list(data[:, 3]/10) + [None]
    xs += list(data[:, 2][[0, -1]])
    ys += list(data[:, 3][[0, -1]]/10)
axes[0].plot(x, y, color="C0", zorder=10)
axes[0].scatter(xs, ys, edgecolor="k", zorder=11, color="C0")
axes[1].plot(x, y, color="C0", zorder=10)
axes[1].scatter(xs, ys, edgecolor="k", zorder=11, color="C0")
skip = 1652
if reactive:
    axes[0].plot(data[:skip, 2], data[:skip, 3]/10, color="C1", zorder=9)
    axes[0].scatter(data[:skip, 2][[0, -1]], data[:skip, 3][[0, -1]]/10, color="C1", zorder=12, edgecolor="k")
    axes[1].plot(data[:skip, 2], data[:skip, 3]/10, color="C1", zorder=9)
    axes[1].scatter(data[:skip, 2][[0, -1]], data[:skip, 3][[0, -1]]/10, color="C1", zorder=12, edgecolor="k")

for i in range(1, 1000):
    path = f"../infretis_{i}.toml"
    if not os.path.exists(path):
        break
    intfs = read_toml(path)["simulation"]["interfaces"]
    c = intfs[-2]

    y0 = c - (1/a)*x0
    axes[0].plot(x0, -y0, ls="--", zorder=3000)
    axes[1].plot(x0, -y0, ls="--", zorder=3000)

toml = "../infretis.toml"
if os.path.exists(toml):
    intfs = read_toml(toml)["simulation"]["interfaces"]
    c = intfs[-2]

    y0 = c - (1/a)*x0
    axes[0].plot(x0, -y0, ls="--", zorder=3000)
    axes[1].plot(x0, -y0, ls="--", zorder=3000)

plt.savefig("progress.png", dpi=150)
plt.show()
