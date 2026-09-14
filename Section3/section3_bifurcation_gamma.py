"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Bifurcation diagram of x_t as a function of gamma (Section 3), holding
beta_x fixed. For each value of gamma in a fine grid, the model is
simulated for T periods and the last T_last periods are plotted against
that gamma value, showing which long-run values of x_t are stable.
Produces Fig. X_bifurcations.

Version 2.4, May 2024
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import math

start = time.time()

def greenGRT_bif5bis(E0, X0, gamma0, betax, betae, alpha1, delta, mu):
    """
    Bifurcation diagram generator: simulates the deterministic model
    (DICE damage function, no gamma/betax switch) for a grid of gamma
    values from gamma_min to gamma_max, and plots the last T_last
    periods of x_t against each gamma value.

    Note: this is a self-contained function distinct from
    section3_baseline_model.greenGR_newnew5bis / section3_scenario_model.
    greenGR_scenario, since it performs a full parameter sweep (a
    separate simulation per gamma value) rather than a single run.

    Parameters
    ----------
    E0, X0 : float
        Initial net emissions and initial relative share of countries acting.
    gamma0 : float
        Initial value of gamma (only used for gamma[0]; overwritten by the
        gamma_range sweep from t=1 onward).
    betax, betae, mu : float
        Peer-pressure, private-utility, and mean-bias parameters (fixed
        across the sweep).
    alpha1 : float
        Effect of participation on emissions growth.
    delta : float
        Sensitivity of gamma to damages (delta=0 recovers the no-damages
        baseline, as used in the published figure).

    Returns
    -------
    E, X, v, gamma, Omega : np.ndarray
        Time series from the LAST gamma value in the sweep only (the
        bifurcation plot itself, built during the loop, is what matters;
        these returned arrays are a by-product, not the main output).
    """

    T=1000
    alpha0=0
    E=np.zeros(T)  #stock greenhouse gasses
    X=np.zeros(T)  #relative share of countries wich take climate action (-1,1)
    v=np.zeros(T)  #value of take action
    VE=np.zeros(T) #growth rate of E
    gamma=np.zeros(T) # intensity of choice 
    Omega=np.zeros(T) #damages function
    sum_E=np.zeros(T) #stock of Emissions
    #for achieving Omega we need
    
    mat=np.zeros(T) # mass of carbon atmosphere 
    Vmat=np.zeros(T) # mass of carbon atmosphere t-1
    
    mup=np.zeros(T)  # mass of upper ocean 
    Vmup=np.zeros(T) # mass of upper ocean t-1
    
    mlo=np.zeros(T)# mass of lower ocean 
    Vmlo=np.zeros(T) # mass of lower ocean t-1
    
    F=np.zeros(T)  #total radiative forcing + FEX exogenous 
    FEX=np.zeros(T)   # exogenous radiative force 
    VFEX=np.zeros(T)
    
    tempAT=np.zeros(T)  #temp atmosphere 
    VtempAT=np.zeros(T) #temp atmosphere t-1
    
    tempLO=np.zeros(T) #temp lower ocean
    VtempLO=np.zeros(T) #temp lower ocean t-1
    
    
    mat[0]=0
    mup[0]=0
    mlo[0]=0
    #F[0]=F0
    tempAT[0]=0
    tempLO[0]=0
    Vmat[0]=3120
    Vmup[0]=5628.8
    Vmlo[0]=36706.7
    FEX[0]=0
    VFEX[0]=0.28
    VtempAT[0]=1
    VtempLO[0]=0.0068
    

    X[0]=X0
    eta=3.8
    xi1=0.027
    xi2=eta/3
    xi3=0.018
    xi4=0.005
    
    phi11=0.9817
    phi21=0.0080
    phi12=0.0183
    phi22=0.9915
    phi23=0.0005
    phi32=0.0001
    phi33=0.9999
    fex=0.005
    MATP=2156.2
    VE[0]=0
    X[0]=X0
    E[0]=E0 #(alpha0-alpha1*X[0])
    gamma[0]=gamma0
    sum_E[0]=1730
    
    gamma_min=0.001
    gamma_max=2
    
    step=0.001
    gamma_range = np.linspace(gamma_min, gamma_max, int(1/step+1))
    Z=np.zeros((int(1/step+1),T))   
    
    flag=-1

    # For each value of gamma in gamma_range, re-simulate the model from
    # scratch and store the full x_t path in Z (one row per gamma value).
    for i in gamma_range:
        flag = flag + 1
        
        
        for t in range(T-1):  
            
            
            #mat[t]=E[t]+phi11*Vmat[t]+phi21*Vmup[t]
            mat[t]=sum_E[t]+phi11*Vmat[t]+phi21*Vmup[t]
            mup[t]=phi12*Vmat[t]+phi22*Vmup[t]+phi32*Vmlo[t]
            mlo[t]=phi23*Vmat[t]+phi33*Vmlo[t] 
            FEX[t]=VFEX[t]+fex
            F[t]=eta*(math.log2(mat[t]/MATP))+FEX[t]  # FEX exogenous ?
            
            tempAT[t]=VtempAT[t]+xi1*(F[t]-xi2*VtempAT[t]-xi3*(VtempAT[t]-VtempLO[t])) #temp atmosphere 
            tempLO[t]=VtempLO[t]+xi4*(VtempAT[t]-VtempLO[t])  #temp lower ocean
            
            Omega[t]=1-1/(1+0.0022*(tempAT[t]**2))  #omega DICE2016R
            
            
            v[t]=gamma[t]*(betax*X[t]+betae*E[t]-mu) 
            E[t+1]=E[t]+E[t]*(alpha0-alpha1*X[t])
            gamma[t+1]=i-delta*Omega[t]     # here is gamma_range
            VE[t+1]=(E[t+1]-E[t])*100/E[t]
            X[t+1]=((np.exp(v[t])-1)/(1+np.exp(v[t])))
            
           
            sum_E[t+1]=sum_E[t]+E[t]
            Vmat[t+1]=mat[t]            
            Vmup[t+1]=mup[t]
            Vmlo[t+1]=mlo[t]
            VtempAT[t+1]=tempAT[t]
            VtempLO[t+1]=tempLO[t]
            VFEX[t+1]=FEX[t]
  
        Z[flag,:] = X

    T_last = 100
    Nr =int( T - T_last)

    # Plot only the last T_last periods of each simulated path: this
    # shows the long-run (converged) values of x_t for each gamma,
    # which is what a bifurcation diagram displays.
    
    plt.style.use('ggplot')
    #plt.style.use('bmh')
    plt.grid(True)
    plt.xlabel(r'$\gamma$')
    plt.xlim(gamma_min, gamma_max)
    #plt.title('Bifurcation diagram')
    #plt.legend(bbox_to_anchor=(0, 1), loc='upper left', borderaxespad=0.1)
    
    for i in range(len(gamma_range)): 
        plt.plot(gamma_range[i]*np.ones((T_last)), Z[i, Nr:T], linestyle='dashdot', linewidth=0.1, markersize=1, marker=',', color='red')
        #plt.title('Bifurcation diagram')
        plt.ylabel(r'$x_t$')
        plt.ylim(-1.1, 1.1)
        plt.xlim(gamma_min, gamma_max)


    
  
    return E, X, v, gamma, Omega

# Bifurcation diagram with the estimated parameters (beta_x=15, delta=0)
# — Fig. X_bifurcations
E, X, v, gamma, Omega=greenGRT_bif5bis(E0=36.8, X0=-.97, gamma0=1, betax=15, betae=0.18, alpha1=0.02, delta=0, mu=9.77)
plt.savefig('bif_gamma_estimated.png', dpi=1200, format='png', bbox_inches='tight')
plt.show()
print(f"Elapsed time: {time.time()-start:.1f} seconds")