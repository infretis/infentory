"""Plot the cosine bump potential in the style of flat_pot.pdf."""
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
})

# Parameters from infretis.toml
bleft  = -0.1
bright =  0.1
nbump  = 1
deltaf = 1.0

# Interfaces: lambda_-1, lambda_0, lambda_1, lambda_2
lambda_m1 = -0.2
intfs = [lambda_m1, -0.1, 0.0, 0.1]

# x grid
x = np.linspace(-0.5, 0.5, 40001)

# Cosine bump potential (no walls: xleft/xright at -100/100)
L      = bright - bleft          # width of one bump: 0.2
bright2 = bleft + nbump * L      # end of bump region: 0.1

pot = np.zeros_like(x)
mask = (x >= bleft) & (x <= bright2)
pot[mask] = deltaf / 2.0 * (1.0 - np.cos(2.0 * np.pi * (x[mask] - bleft) / L))

# --- Figure ---
fig, ax = plt.subplots(figsize=(3.5, 2.8))

xlim = [-0.45, 0.45]
ylim = [-0.15, 1.35]

# Grey shading: A spans xlim[0] to lambda_0 (lambda_-1 is inside A)
ax.axvspan(xlim[0], intfs[1],  facecolor='grey', alpha=0.2)
ax.axvspan(intfs[-1], xlim[1], facecolor='grey', alpha=0.2)

# Interface labels — placed in a gap in the vlines
y_label   = 1.20          # vertical centre of label text
y_gap_lo  = 1.13          # vline stops here below the label
y_gap_hi  = 1.28          # vline resumes here above the label

# Dashed vlines with a gap for the text
ax.vlines(intfs, ylim[0], y_gap_lo, colors='black', linestyle='--', linewidth=0.6)
ax.vlines(intfs, y_gap_hi, ylim[1], colors='black', linestyle='--', linewidth=0.6)

# Potential curve
ax.plot(x, pot, linewidth=1.5, color='C0')

# A / B labels (A is centred in the full grey band from xlim[0] to lambda_0)
ax.text((xlim[0] + intfs[1]) / 2+0.03, ylim[1] - 0.65, '$A$',
        ha='center', va='center', fontsize=16)
ax.text(intfs[-1] + 0.09, ylim[1] - 0.65, '$B$',
        ha='center', va='center', fontsize=16)

# Interface labels inside the vline gap
intf_labels = [r'$\lambda_{-1}$', r'$\lambda_0$', r'$\lambda_1$', r'$\lambda_2$']
for intf, label in zip(intfs, intf_labels):
    ax.text(intf, y_label, label, ha='center', va='center', fontsize=11)

ax.set_xlim([-0.3, 0.3])
ax.set_ylim(ylim)
ax.minorticks_off()
ax.set_xticks(np.arange(-0.3, 0.31, 0.2))
ax.set_xlabel(r'Position $x$ [red.]', fontsize=12)
ax.set_ylabel(r'Potential $V(x)$ [$k_B T$]', fontsize=12)
ax.set_title('cosine bump', fontsize=15)

plt.tight_layout()
plt.savefig('cos_bump_pot.pdf', dpi=300, bbox_inches='tight')
plt.savefig('cos_bump_pot.png', dpi=300, bbox_inches='tight')
print("Saved cos_bump_pot.pdf and cos_bump_pot.png")
