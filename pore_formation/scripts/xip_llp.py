import matplotlib.pyplot as plt
import numpy as np
import os
from inftools.misc.infinit_helper import read_toml

a = 0.25

x0 = np.linspace(0, 3, 100)
c = [-1.2, 3.5]
for c0 in c:
    y0 = c0 - (1/a)*x0
    plt.plot(x0, y0, ls="--", zorder=3, color="k")

x2 = np.loadtxt("y.txt")
y2 = np.loadtxt("x.txt")/10
ufree2 = np.loadtxt("unifree.txt")
plt.contour(x2, -y2, ufree2 - np.min(ufree2), colors='black', alpha=0.25, zorder=1, levels=20)

x, y, xs, ys = [], [], [], []
for i in range(0, 1000):
    path = f"../load/{i}/order.txt"
    if not os.path.exists(path):
        break
    data = np.loadtxt(path)
    x += list(data[:, 2]) + [None]
    y += list(-data[:, 3]/10) + [None]
    xs += list(data[:, 2][[0, -1]])
    ys += list(-data[:, 3][[0, -1]]/10)
plt.plot(x, y)
plt.scatter(xs, ys, edgecolor="k", zorder=4)

for i in range(1, 1000):
    path = f"../infretis_{i}.toml"
    if not os.path.exists(path):
        break
    intfs = read_toml(path)["simulation"]["interfaces"]
    c = intfs[-2]

    y0 = c - (1/a)*x0
    plt.plot(x0, y0, ls="--", zorder=3)

toml = "../infretis.toml"
if os.path.exists(toml):
    intfs = read_toml(toml)["simulation"]["interfaces"]
    c = intfs[-2]

    y0 = c - (1/a)*x0
    plt.plot(x0, y0, ls="--", zorder=3)

plt.ylim([-3.5, -1])
plt.xlim([0.2, 2.2])

plt.ylabel(r"- $d_{\rm LLP}$ [nm]")
plt.xlabel(r"$\xi_{\rm p}$")

plt.show()
