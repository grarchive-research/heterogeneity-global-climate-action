"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Phase plots of x_t (Section 3, Appendix): x_t plotted against x_{t-1}
for different values of gamma (Fig. Phaseplot1) and beta_x (Fig.
Phaseplot2), holding the other baseline parameters fixed. The gamma=1,
beta_x=1.1 panel (Phase_plot_gamma_1.png) is the baseline case and is
shared between the two figures in the paper.

Version 2.4, May 2024
"""

import matplotlib.pyplot as plt
from section3_baseline_model import greenGR_newnew5bis


def plot_phase(gamma0, betax, filename):
    """
    Simulate the baseline model with the given gamma0/betax and save a
    phase plot of x_t against x_{t-1} (periods 29 to 219, discarding the
    initial transient).

    Parameters
    ----------
    gamma0 : float
        Degree of heterogeneity used for this panel.
    betax : float
        Peer-pressure parameter used for this panel.
    filename : str
        Output PNG filename (relative path, no subfolder).
    """
    E, X, v, gamma, Omega = greenGR_newnew5bis(
        E0=36.8, X0=-.97, gamma0=gamma0, betax=betax, betae=0.18,
        alpha1=0.02, delta=0, mu=9.77, a=0
    )

    X_lagged = X[29:220]    # x_t
    X_original = X[28:219]  # x_{t-1}

    plt.figure()
    plt.plot(X_original, X_lagged)
    plt.xlabel(r"$x_{t-1}$")
    plt.ylabel(r"$x_{t}$")
    plt.savefig(filename, dpi=1200, format='png', bbox_inches='tight')
    plt.show()


if __name__ == "__main__":

    # --- Fig. Phaseplot1: phase plot for different values of gamma (beta_x=1.1 fixed) ---
    plot_phase(gamma0=1,   betax=1.1, filename='Phase_plot_gamma_1.png')    # baseline, reused in Phaseplot2
    plot_phase(gamma0=1.5, betax=1.1, filename='Phase_plot_gamma_1_5.png')
    plot_phase(gamma0=2,   betax=1.1, filename='Phase_plot_gamma_2.png')

    # --- Fig. Phaseplot2: phase plot for different values of beta_x (gamma=1 fixed) ---
    # Phase_plot_gamma_1.png above already covers beta_x=1.1
    plot_phase(gamma0=1, betax=3,  filename='Phase_plot_beta_3.png')
    plot_phase(gamma0=1, betax=15, filename='Phase_plot_beta_15.png')
