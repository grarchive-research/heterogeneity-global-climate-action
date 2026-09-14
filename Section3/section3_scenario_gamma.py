"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Gamma scenario (Section 3): compares the baseline (gamma=1, no switch)
against two alternative degrees of heterogeneity applied after
t=switch_t (2019). Produces Fig. XandE_gamma — illustrative of the
model's qualitative dynamics, not a forecast of actual future outcomes.

Version 2.4, May 2024
"""

import matplotlib.pyplot as plt
from section3_scenario_model import greenGR_scenario

# --- Panel 1: participation x_t under three gamma scenarios (Fig. XandE_gamma) ---


# Scenario: gamma switches from 1 (baseline) to 3 after 2019
E, X, v, gamma, Omega= greenGR_scenario(E0=36.8, X0=-.97, gamma0=1, betax=1.1, gamma2=3, betae=0.18, alpha1=0.02, delta=0, mu=9.77, a=0)
plt.grid(True)
plt.xlabel(r'$t$')
plt.ylim([-1.1, 1.1])
plt.plot(X, color = 'blue', linewidth=1.0, linestyle='dotted', label=r'$x_t$')
plt.xlim([30, 60])
plt.xticks([0,10,20,30,40,50,60], ['1989','1999','2009','2019','2029','2039','2049'])
plt.style.use('ggplot')

# Scenario: gamma switches from 1 (baseline) to 15 after 2019
E, X, v, gamma, Omega= greenGR_scenario(E0=36.8, X0=-.97, gamma0=1, betax=1.1, gamma2=15, betae=0.18, alpha1=0.02, delta=0, mu=9.77, a=0)
plt.grid(True)
plt.xlabel(r'$t$')
plt.ylim([-1.1, 1.1])
plt.plot(X, color = 'black', linewidth=1.0, linestyle='--', label=r'$x_t$')
plt.xlim([30, 60])
plt.xticks([0,10,20,30,40,50,60], ['1989','1999','2009','2019','2029','2039','2049'])
plt.style.use('ggplot')

# Baseline scenario: gamma stays at 1 throughout (no switch)
E, X, v, gamma, Omega= greenGR_scenario(E0=36.8, X0=-.97, gamma0=1, betax=1.1, gamma2=1, betae=0.18, alpha1=0.02, delta=0, mu=9.77, a=0)
plt.grid(True)
color = 'tab:red'
plt.xlabel(r'$t$')
plt.ylabel(r'$x_t$', color = color)
plt.ylim([-1.1, 1.1])
plt.plot(X, color = color, linewidth=1.0, linestyle='dashed', label=r'$x_t$')
plt.xlim([0, 60])
plt.xticks([0,10,20,30,40,50,60], ['1989','1999','2009','2019','2029','2039','2049'])
plt.tick_params(axis ='y', labelcolor = color)
plt.style.use('ggplot')
plt.savefig('X_gamma_scenario.png', dpi=1200, format='png', bbox_inches='tight')
plt.show()

# --- Panel 2: net emissions E_t under the same three gamma scenarios ---

# Baseline scenario: gamma stays at 1 throughout (no switch — baseline reference)
E, X, v, gamma, Omega= greenGR_scenario(E0=36.8, X0=-.97, gamma0=1, betax=1.1, gamma2=1, betae=0.18, alpha1=0.02, delta=0, mu=9.77, a=0)
plt.grid(True)
color = 'tab:green'
plt.xlabel(r'$t$')
plt.ylabel(r'$E_t$', color = color)
plt.ylim([0, 100])
plt.plot(E, color = color, linewidth=1.0, linestyle='dashed', label=r'$E_t$')
plt.xlim([0, 60])
plt.xticks([0,10,20,30,40,50,60], ['1989','1999','2009','2019','2029','2039','2049'])
plt.tick_params(axis ='y', labelcolor = color)
plt.style.use('ggplot')

# Scenario: gamma switches from 1 to 3 after 2019
E, X, v, gamma, Omega= greenGR_scenario(E0=36.8, X0=-.97, gamma0=1, betax=1.1, gamma2=3, betae=0.18, alpha1=0.02, delta=0, mu=9.77, a=0)
plt.grid(True)
plt.xlabel(r'$t$')
plt.plot(E, color = 'blue', linewidth=1.0, linestyle='dotted', label=r'$E_t$')
plt.xlim([30, 60])
plt.xticks([0,10,20,30,40,50,60], ['1989','1999','2009','2019','2029','2039','2049'])
plt.style.use('ggplot')

# Scenario: gamma switches from 1 to 15 after 2019
E, X, v, gamma, Omega= greenGR_scenario(E0=36.8, X0=-.97, gamma0=1, betax=1.1, gamma2=15, betae=0.18, alpha1=0.02, delta=0, mu=9.77, a=0)
plt.grid(True)
plt.xlabel(r'$t$')
plt.plot(E, color = 'black', linewidth=1.0, linestyle='--', label=r'$E_t$')
plt.xlim([30, 60])
plt.xticks([0,10,20,30,40,50,60], ['1989','1999','2009','2019','2029','2039','2049'])
plt.style.use('ggplot')
plt.savefig('E_gamma_scenario.png', dpi=1200, format='png', bbox_inches='tight')
plt.show()


