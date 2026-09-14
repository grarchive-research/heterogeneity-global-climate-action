"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Compares the deterministic baseline model (see section3_baseline_model.py)
against observed historical data on climate participation (x_t) and net
GHG emissions (E_t), 1989-2021. Produces the model-data comparison figure
for Section 3 of "Heterogeneity and Global Climate Action"
(Fig. XandE_actual_simulated).

Data source: xt_GHG_series.xlsx (Climate Laws of the World Database;
Our World in Data)

Version 2.4, May 2024
"""


import numpy as np 
import matplotlib.pyplot as plt
import math
import pandas as pd
from section3_baseline_model import greenGR_newnew5bis

# NOTE: 'year' is currently unused in this script (kept for reference/future use)
year=np.arange(1989,2989,1)

# Simulate the baseline model with the estimated parameters (Table 2, appendix)
E, X, v, gamma, Omega= greenGR_newnew5bis(E0=36.8, X0=-.97, gamma0=1, betax=1.1, betae=0.18, alpha1=0.02, delta=0, mu=9.77, a=0)

# Load observed historical data: x_t (participation) and GHG (net emissions)
file_path = 'xt_GHG_series.xlsx'
df = pd.read_excel(file_path)
print(df)

# Quick sanity-check scatter: relationship between observed participation and emissions
df.plot(x='x_t', y='GHG', kind='scatter')

plt.title('Scatter Plot of $x_t$ vs $GHG$')
plt.xlabel('x_t')
plt.ylabel('GHG')
plt.show()

# Observed time series: participation (x_t) and emissions (GHG) over 1989-2021
df.plot(x='date', y='GHG', kind='line', color = 'green')
df.plot(x='date', y='x_t', kind='line', color= 'red')
plt.xlabel('Time')
plt.ylabel('Value')
plt.show()

# --- Participation: observed vs simulated (Fig. XandE_actual_simulated, panel a) ---
df['participation'] = X[0:len(df)]
print(df)

plt.plot(df['date'], df['x_t'], label='Actual', color='blue')
plt.plot(df['date'], df['participation'], label='Similated', color='red')
plt.xlabel('$t$')
plt.ylim([-1.1, 1.1])
plt.ylabel(r'Participation')
plt.savefig('X_actual_simulated.png', dpi=1200, format='png', bbox_inches='tight')

# --- Net emissions: observed vs simulated (Fig. XandE_actual_simulated, panel b) ---
df['E_t'] = E[0:len(df)]
print(df)


plt.plot(df['date'], df['GHG'], label='Actual', color='blue')
plt.plot(df['date'], df['E_t'], label='Simulated', color='green')

plt.xlabel('Time')
plt.ylim([0, 100])
plt.ylabel('Net Emissions')
plt.savefig('E_actual_simulated.png', dpi=1200, format='png', bbox_inches='tight')


