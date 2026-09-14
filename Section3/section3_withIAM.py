"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Compares the deterministic baseline model (see section3_baseline_model.py)
under three damage-function assumptions — no damages (baseline, delta=0),
DICE2016R (a=0, delta=0.99) and Weitzman (2012) (a=1, delta=0.99) — for
Section 3 of "Heterogeneity and Global Climate Action". Produces
Fig. Comparison_evolution_all, panel a (participation x_t, estimated
parameters, gamma0=1).

Version 2.4, May 2024
"""

import matplotlib.pyplot as plt
from section3_baseline_model import greenGR_newnew5bis as greenGR_newnew6witIAM
 

# --- Fig. Comparison_evolution_all, panel a: x_t under baseline / DICE / Weitzman ---

# DICE damage function (a=0), endogenous gamma via damages (delta=0.99)
E, X, v, gamma, Omega =greenGR_newnew6witIAM(E0=36.8, X0=-.97, gamma0=1, betax=1.1, betae=0.18, alpha1=0.02, delta=0.99, mu=9.77, a=0)
plt.ylabel('$x_t$')
plt.xlabel('$t$')
plt.ylim([-1.1, 1.1])
plt.style.use('ggplot')
plt.grid(True)
plt.xlim([0,60])
plt.xticks([0,10,20,30,40,50,60], ['1989','1999','2009','2019','2029','2039','2049'])
plt.plot(X, color ='blue', linewidth=1.0, linestyle='dashed', label=r'$DICE$')

# Weitzman (2012) damage function (a=1), endogenous gamma via damages (delta=0.99)
E, X, v, gamma, Omega =greenGR_newnew6witIAM(E0=36.8, X0=-.97, gamma0=1, betax=1.1, betae=0.18, alpha1=0.02, delta=0.99, mu=9.77, a=1)
plt.plot(X, color ='red', linewidth=1.0, linestyle='dashed', label=r'$Weitzman$')

# No-damages baseline (delta=0), gamma exogenous — same as section3_baseline_model.py
E, X, v, gamma, Omega =greenGR_newnew6witIAM(E0=36.8, X0=-.97, gamma0=1, betax=1.1, betae=0.18, alpha1=0.02, delta=0, mu=9.77, a=0)
plt.plot(X, color ='black', linewidth=1.0, linestyle='dashed', label=r'$Baseline$')

plt.savefig('Comparison_evolution_all_a.png', dpi=1200, format='png', bbox_inches='tight')
plt.show()




