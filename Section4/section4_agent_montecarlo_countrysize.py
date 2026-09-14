"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Monte Carlo pipeline for the agent-based model (Section 4), GDP-share
weights: 2N=194 countries on a fully-connected network, weighted by
relative GDP share (weights_array from dataset.py). Nested double
Monte Carlo — 100 draws of the idiosyncratic error term (inner loop) x
194 starting countries, one per possible first mover (outer loop) —
averaged to produce the mean and standard-deviation band of
participation and net emissions. Produces the published Fig.
Comparison_evolution_discrete_size.

Version 3.0, Jan-Feb 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from dataset import dataset
from section4_agent_baseline import GRT_AMBnetEJ

weights_array, GHG_array = dataset()

    
def MonteCarlo1(N, T, MC_runs):

    E_results = np.zeros((N, MC_runs, T))  # Stores emissions for each country and MC run
    X_results = np.zeros((N, MC_runs, T))  # Stores action shares
    VE_results = np.zeros((N, MC_runs, T))  # Stores emission growth rates
    nC_results = np.zeros((N, MC_runs, T))  # Stores number of acting countries

    for k in range(N):  # Iterate over each country as the starting point
        print(f"Running model with initial country {k}")

        for mc in range(MC_runs):  # Monte Carlo loop for error term variations
            np.random.seed(mc)  # set the seed
            E, X, pi, VE, nC, y, Z, ZZ, A, C, D = GRT_AMBnetEJ(T=T, betan=1.1, betae=0.18, betad=0, N=N, a=0, mu=9.77, start_idx=k, weights_array=weights_array)

            E_results[k, mc, :] = E
            X_results[k, mc, :] = X
            VE_results[k, mc, :] = VE
            nC_results[k, mc, :] = nC

    # Compute mean and standard deviation over all MC runs and starting countries
    mean_E = np.mean(E_results, axis=(0,1))
    mean_X = np.mean(X_results, axis=(0,1))
    mean_VE = np.mean(VE_results, axis=(0,1))
    mean_nC = np.mean(nC_results, axis=(0,1))
    std_E = np.std(E_results, axis=(0,1))
    std_X = np.std(X_results, axis=(0,1))
    std_VE = np.std(VE_results, axis=(0,1))
    std_nC = np.std(nC_results, axis=(0,1))

    return mean_E, mean_VE, mean_X, mean_nC, std_E, std_X, std_VE, std_nC

mean_E, mean_VE, mean_X, mean_nC, std_E, std_X, std_VE, std_nC = MonteCarlo1(194, 500, MC_runs=50)


fig, ax1 = plt.subplots()
T=500
# Plot mean X_t with standard deviation shading
color = 'tab:red'
ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$x_t$', color=color)
ax1.plot(mean_X, color=color, linewidth=1.0, linestyle='dashed', label=r'Mean $x_t$')
ax1.fill_between(range(T), mean_X - std_X, mean_X + std_X, color=color, alpha=0.2)  # Light shading
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim([-1, 1])  # Set limits for x_t
ax2 = ax1.twinx()

# Plot mean E_t with standard deviation shading
color = 'tab:green'
ax2.set_ylabel(r'$E_t$', color=color)
ax2.plot(mean_E, color=color, linewidth=1.0, linestyle='dashed', label=r'Mean $E_t$')
ax2.fill_between(range(T), mean_E - std_E, mean_E + std_E, color=color, alpha=0.2)  # Light shading
ax2.tick_params(axis='y', labelcolor=color)

# Styling
plt.style.use('ggplot')
plt.grid()
plt.xlim([0, 60]) 
plt.xticks([0,10,20,30,40,50,60], ['1989','1999','2009','2019','2029','2039','2049'])
plt.ylim([20, 80])
fig.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2, borderaxespad=0.21)
plt.savefig('MonteCarlo_withcountriesweights.png', dpi=1200, format='png', bbox_inches='tight')
plt.show()




