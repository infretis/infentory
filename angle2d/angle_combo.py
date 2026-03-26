import numpy as np
import os
from inftools.tistools.max_op import COLS
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import scienceplots
# apply science style yeah
plt.style.use('science')

# set colors
start_color = "#ff7f0e"  # blue
end_color   = "#1f77b4"  # orange
cmap = LinearSegmentedColormap.from_list(
    "blue_to_orange",
    [end_color, start_color]
)
colors = cmap(np.linspace(1, 0, len(np.arange(0, 91, 5))))[::-1]

# define figure
fig, axs = plt.subplots(2, 3, figsize=(8.5, 4.0), layout="constrained")

# set angle and collect conditional and unconditional
angles = np.linspace(0, 90, 19)
ufreesb = []
cfreef = []
cfreeb = []

# comparison points
cmp_f =  0.7
cmp_b = -0.7

# transitions state through the x axis, it is not zero btw
x_ts = -0.06275

# Read in datas
ufreexes = []
ufreeyes = []
for angle in angles:
    ufreex = np.loadtxt(f"run_f1/wham_00/histo_xval.txt")
    ufreey = np.loadtxt(f"frees/ufree_{angle:02.0f}.txt")
    ufreexes.append(ufreex)
    ufreeyes.append(ufreey)

    index = np.argmin(np.abs(ufreex-x_ts))
    ufreesb.append(ufreey[index])

    cfreefx = np.loadtxt(f"run_f1/wham_{angle:02.0f}/histo_xval.txt")
    cfreefy = np.loadtxt(f"run_f1/wham_{angle:02.0f}/histo_free_energy.txt")
    index = np.argmin(np.abs(cfreefx-cmp_f))
    cfreef.append(cfreefy[index])

    cfreebx = np.loadtxt(f"run_b1/wham_{angle:02.0f}/histo_xval.txt")
    cfreeby = np.loadtxt(f"run_b1/wham_{angle:02.0f}/histo_free_energy.txt")
    index = np.argmin(np.abs(cfreebx-cmp_b))
    cfreeb.append(cfreeby[index])

#  plot 2d conditional free energies
for idx, path, vmax in zip([0, 2], ["run_f1", "run_b1"], [24, 33]):
    histx = np.loadtxt(f"{path}/wham2d/histo_xval.txt")
    histy = np.loadtxt(f"{path}/wham2d/histo_xval.txt")
    histf = np.loadtxt(f"{path}/wham2d/histo_free_energy.txt")
    axs[0][idx].pcolormesh(histx, histy, histf.T, rasterized=True, vmin=0, vmax=vmax,)
    axs[0][idx].set_xlim([-1.8, 1.8])
    axs[0][idx].set_ylim([-1.8, 1.8])
    contour = axs[0][idx].contour(histx, histy, histf.T, levels=15, colors='k', linewidths=0.5, alpha=0.1, vmin=-2, vmax=vmax, zorder=100)

# plot the potential grid
def pot(x, y, a=1.3128, b=2, c=0.25):
    return a*(y**4 + x**4 -b*x**2 - c*x)

## Create grid
x = np.linspace(-2.2, 2.2, 200)
y = np.linspace(-2.2, 2.2, 200)
X, Y = np.meshgrid(x, y)
Z = pot(X, Y)
contour = axs[0][0].contour(X, Y, Z, levels=25, colors='k', linewidths=0.5, alpha=0.5, vmin=0, vmax=30, zorder=100)
contour = axs[0][2].contour(X, Y, Z, levels=25, colors='k', linewidths=0.5, alpha=0.5, vmin=0, vmax=30, zorder=100)

# plot the delta free energies
skip0 = 11
axs[0][1].plot(angles, ufreesb - np.min(ufreeyes[0]), color="k", lw=1.5, zorder=1)
axs[0][1].plot(angles[:skip0], ufreesb[:skip0] - np.min(ufreeyes[0][ufreexes[0]<0]),  color="k", lw=1.5, zorder=1)
axs[0][1].plot(angles, cfreef-np.min(cfreefy), color=COLS[0], lw=1.5, zorder=2)
axs[0][1].plot(angles, cfreeb-np.min(cfreeby), color=COLS[3], lw=1.5, zorder=3)
axs[0][1].scatter(angles[:skip0], ufreesb[:skip0]-np.min(ufreeyes[0][ufreexes[0]<0]), color="k", s=20, zorder=1, marker="o", label=r"$F_{\rightarrow}$")

axs[0][1].scatter(angles, ufreesb - np.min(ufreeyes[0]), color="k", s=20, zorder=1, marker="x",label=r"$F_{\leftarrow}$")
axs[0][1].scatter(angles, cfreef-np.min(cfreefy), color=COLS[0], s=20, zorder=2, marker="o",label=r"$F_{\mathcal{A}\rightarrow}$")
axs[0][1].scatter(angles, cfreeb-np.min(cfreeby), color=COLS[3], s=20, zorder=3, marker="x",label=r"$F_{\mathcal{B}\leftarrow}$")
axs[0][1].legend(handlelength=0.2, loc="lower left", ncol=2, columnspacing=0.85, handletextpad=0.5, fontsize=9)

axs[0][0].set_xlabel("x")
axs[0][0].set_ylabel("y")
axs[0][2].set_xlabel("x")
axs[0][2].set_ylabel("y")
axs[0][1].set_xlabel(r"$\theta$")
axs[0][1].set_ylabel(r"$\Delta F(\theta)$ $(k_BT)$")

# plot the inset color bar in [1][0]
norm2 = mcolors.Normalize(vmin=0, vmax=90)
sm2 = cm.ScalarMappable(norm=norm2, cmap=cmap)
sm2.set_array([])

cax = fig.add_axes([0.08, 0.170, 0.05, 0.010])
cbar = plt.colorbar(sm2, cax=cax, orientation="horizontal")
cbar.set_ticks([0, 90])
cbar.set_ticklabels([r"$0^\circ$", r"$90^\circ$"])
cbar.set_label(r"$\theta$", labelpad=-8)

simkb = 0.07
# the outliers of free energy is not smooth, just skip them
skip = 15
# The energy difference between the two minima across X
delta = -1.6459725407647032 + 0.9898830229008313

# plot the conditional and unconditional 1D free energies
for idx, angle in enumerate(angles):
    datax = ufreexes[idx]
    dataf = ufreeyes[idx]
    alpha = 0.2 if idx not in (0, len(angles)-1) else 1.0
    axs[1][1].plot(datax[dataf!=np.inf][skip:-skip], dataf[dataf!=np.inf][skip:-skip], zorder=3, alpha=alpha, color=colors[idx], lw=1.5)
    dataff = np.loadtxt(f"run_f1/wham_{angle:02.0f}/histo_free_energy.txt")
    axs[1][0].plot(datax[dataff!=np.inf][skip:-skip], dataff[dataff!=np.inf][skip:-skip] - delta/simkb, zorder=3, alpha=alpha, color=colors[idx], lw=1.5)
    datafb = np.loadtxt(f"run_b1/wham_{angle:02.0f}/histo_free_energy.txt")
    axs[1][2].plot(datax[datafb!=np.inf][skip:-skip], datafb[datafb!=np.inf][skip:-skip], zorder=3, alpha=alpha, color=colors[idx], lw=1.5)

# mark the comparison points
axs[1][0].axvline(cmp_f, color="r", ls="--", lw="1.5",    alpha=0.2)
axs[1][1].axvline(x_ts, color="r", ls="--", lw="1.5", alpha=0.2)
axs[1][2].axvline(cmp_b, color="r", ls="--", lw="1.5",    alpha=0.2)

axs[0][1].set_xticks([0, 15, 30, 45, 60, 75, 90])
axs[0][1].set_yticks([0,10, 20, 30])

# plot the actual potential 1d slices
datax = np.linspace(-2, 2, 100)
free_y = pot(0, datax)
free_x = pot(datax, 0)
for i in range(3):
    delta2 = delta if i == 0 else 0
    axs[1][i].plot(datax[14:], (free_x[14:] - np.min(free_x[np.sum(datax>0):]))/simkb, color="k", zorder=1, lw=1.5, ls="-")
    axs[1][i].plot(datax, (free_y - delta2)/simkb, color="k", zorder=1, lw=1.5, ls="--")
    axs[1][i].set_xlim([-1.8, 1.8])
    axs[1][i].set_ylim([-2.5, 32])

    axs[1][i].set_xlabel(r"$\chi(x, y, \theta) = x\cos\theta+y\sin\theta$")
    axs[1][i].set_ylabel(r"Energy ($k_BT$)")


# draw the theta..
from matplotlib.patches import Arc, FancyArrowPatch, FancyBboxPatch
theta = 60
x0=0.05
y0=0.05
L = 0.175
dx = L * np.cos(np.deg2rad(theta))
dy = L * np.sin(np.deg2rad(theta))
bg = FancyBboxPatch(
    (x0 - 0.01, y0 - 0.01),
    width=0.15,
    height=0.15,
    boxstyle="round,pad=0.02,rounding_size=0.03",
    transform=axs[0][0].transAxes,
    facecolor="white",
    edgecolor="none",
    zorder=200,
    alpha=0.8,
)
axs[0][0].add_patch(bg)
arc = Arc(
    (x0, y0),
    width=0.125,
    height=0.15,
    angle=0,
    theta1=8,
    theta2=theta,
    color="k",
    linewidth=0.75,
    transform=axs[0][0].transAxes,
    zorder=200,
)
axs[0][0].add_patch(arc)
axs[0][0].text(
    -1.25,
    -1.40,
    r"$\theta$",
    color="k",
    fontsize=10,
    ha="center", va="center",
    zorder=200,
)
arrow = FancyArrowPatch(
    (x0, y0), (x0 + dx, y0 + dy),
    arrowstyle="->",
    color="k",
    linewidth=0.75,
    mutation_scale=12,
    transform=axs[0][0].transAxes,
    zorder=200,
)
axs[0][0].add_patch(arrow)
arrow = FancyArrowPatch(
    (x0-0.005, y0+0.01), (x0 + dx*1.2, y0+0.01),
    arrowstyle="-",
    color="k",
    linewidth=0.75,
    mutation_scale=12,
    transform=axs[0][0].transAxes,
    zorder=200,
)
axs[0][0].add_patch(arrow)

axs[1][0].text(
    0.550, 0.885, r"$F_\mathcal{A}(x)$",
    transform=axs[1][0].transAxes,
    ha="center", va="center", fontsize=9,
)
axs[1][0].text(
    0.500, 0.250, r"$F_\mathcal{A}(y)$",
    transform=axs[1][0].transAxes,
    ha="center", va="center", fontsize=9,
)

axs[1][1].text(
    0.385, 0.80, r"$F(x)$",
    transform=axs[1][1].transAxes,
    ha="center", va="center", fontsize=9,
)
axs[1][1].text(
    0.220, 0.150, r"$F(y)$",
    transform=axs[1][1].transAxes,
    ha="center", va="center", fontsize=9,
)

axs[1][2].text(
    0.45, 0.885, r"$F_\mathcal{B}(x)$",
    transform=axs[1][2].transAxes,
    ha="center", va="center", fontsize=9,
)
axs[1][2].text(
    0.220, 0.150, r"$F_\mathcal{B}(y)$",
    transform=axs[1][2].transAxes,
    ha="center", va="center", fontsize=9,
)

axs[0][0].set_title(r"$F_\mathcal{A}(x, y)$", x=0.5, y=0.775, bbox=dict(facecolor="white", edgecolor="none", boxstyle="round,pad=0.25", alpha=0.75), zorder=400, fontsize=9)
axs[0][2].set_title(r"$F_\mathcal{B}(x, y)$",  x=0.5, y=0.775, bbox=dict(facecolor="white", edgecolor="none", boxstyle="round,pad=0.25", alpha=0.75), zorder=400, fontsize=9)


from matplotlib.transforms import ScaledTranslation
labels = [r"(a)", r"(b)", r"(c)", r"(d)", r"(e)", r"(f)"]

for idx, (ax, label) in enumerate(zip(axs.flatten(), labels)):
    ax.text(
        0.0, 1.0, label, transform=(
            ax.transAxes + ScaledTranslation(-32/72, -12/72, fig.dpi_scale_trans)),
        fontsize='medium', va='bottom', fontfamily='serif',)

# plt.savefig("angle_combo.pdf", bbox_inches="tight", dpi=300)
plt.show()
