"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Robustness check figures for autonomous emissions growth (alpha0), for
the supplementary appendix of "Heterogeneity and Global Climate Action".
Compares the alpha0=0 baseline (section3_baseline_model.py) against the
point estimate alpha0=0.008, alpha1=0.014 (Table 1, column 1, with
constant), and maps long-run participation x* over a grid of alpha0/alpha1
values within two standard errors of the estimates.

Version 1.0, September 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from section3_baseline_model import greenGR_newnew5bis as greenGR_baseline
from section3_baseline_model_alpha import greenGR_newnew5bis as greenGR_alpha0

plt.style.use('ggplot')


# ============ Figure A: x_t and E_t, baseline vs autonomous-growth case ============

# Baseline: alpha0=0, alpha1=0.02 (Table 1, column 2, no constant)
E_base, X_base, v, gamma, Omega = greenGR_baseline(E0=36.8, X0=-.97, gamma0=1, betax=1.1, betae=0.18,
                                                     alpha1=0.02, delta=0, mu=9.77, a=0)

# Autonomous growth: alpha0=0.008, alpha1=0.014 (Table 1, column 1, with constant)
E_new, X_new, v2, gamma2, Omega2 = greenGR_alpha0(E0=36.8, X0=-.97, gamma0=1, betax=1.1, betae=0.18,
                                                    alpha1=0.014, delta=0, mu=9.77, a=0, alpha0=0.008)

fig = plt.figure(figsize=(10, 4.2))

ax1 = fig.add_subplot(1, 2, 1)
ax1.plot(X_base, color='tab:red', linewidth=1.2, linestyle='solid', label=r'Baseline ($\alpha_0=0$)')
ax1.plot(X_new, color='tab:red', linewidth=1.2, linestyle='dashed', label=r'$\alpha_0=0.008$')
ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$x_t$')
ax1.set_ylim([-1.1, 1.1])
ax1.set_xlim([0, 60])
ax1.set_xticks([0, 10, 20, 30, 40, 50, 60])
ax1.set_xticklabels(['1989', '1999', '2009', '2019', '2029', '2039', '2049'])
ax1.grid(True)
ax1.legend(loc='lower right', fontsize=8)

ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(E_base, color='tab:green', linewidth=1.2, linestyle='solid', label=r'Baseline ($\alpha_0=0$)')
ax2.plot(E_new, color='tab:green', linewidth=1.2, linestyle='dashed', label=r'$\alpha_0=0.008$')
ax2.set_xlabel(r'$t$')
ax2.set_ylabel(r'$E_t$')
ax2.set_xlim([0, 60])
ax2.set_xticks([0, 10, 20, 30, 40, 50, 60])
ax2.set_xticklabels(['1989', '1999', '2009', '2019', '2029', '2039', '2049'])
ax2.grid(True)
ax2.legend(loc='lower right', fontsize=8)

plt.tight_layout()
plt.savefig('X_E_alpha0_robustness.png', dpi=1200, format='png', bbox_inches='tight')
plt.show()

# ============ Figure B: sensitivity heatmap over the +-2SE range ============
# alpha0 = 0.008 (se=0.005), alpha1 = 0.014 (se=0.006) -- Table 1, column 1
alpha0_range = np.linspace(0.008 - 2*0.005, 0.008 + 2*0.005, 9)
alpha1_range = np.linspace(0.014 - 2*0.006, 0.014 + 2*0.006, 9)

Xstar_grid = np.full((len(alpha0_range), len(alpha1_range)), np.nan)

for i, a0 in enumerate(alpha0_range):
    for j, a1 in enumerate(alpha1_range):
        if a1 <= 0:
            continue
        E, X, v, gamma, Omega = greenGR_alpha0(E0=36.8, X0=-.97, gamma0=1, betax=1.1, betae=0.18,
                                                alpha1=a1, delta=0, mu=9.77, a=0, alpha0=a0, T=2000)
        if np.isfinite(X[-1]) and np.isfinite(E[-1]) and E[-1] < 200:
            Xstar_grid[i, j] = X[-1]
        # else: leave as NaN -> divergent / no stable equilibrium

fig2, ax = plt.subplots(figsize=(6, 5))
cmap = plt.cm.RdBu_r
cmap.set_bad(color='lightgray')
im = ax.imshow(Xstar_grid, origin='lower', aspect='auto', cmap=cmap, vmin=-1, vmax=1,
               extent=[alpha1_range[0], alpha1_range[-1], alpha0_range[0], alpha0_range[-1]])
ax.set_xlabel(r'$\alpha_1$')
ax.set_ylabel(r'$\alpha_0$')
cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r'long-run $x^*$')
ax.plot(0.014, 0.008, marker='*', color='black', markersize=14)
ax.annotate('point estimate', (0.014, 0.008), textcoords="offset points", xytext=(8, 8), fontsize=8)
plt.tight_layout()
plt.savefig('alpha0_sensitivity_heatmap.png', dpi=1200, format='png', bbox_inches='tight')
plt.show()