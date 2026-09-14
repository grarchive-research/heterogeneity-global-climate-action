"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Deterministic aggregate model (Section 3) with a one-time parameter
switch at t=switch_t. Generalizes section3_baseline_model.py to support
scenarios where beta_x or gamma change after a given period (e.g. after
2019). Used to produce the beta_x and gamma conditional baseline
scenarios (Fig. XandE_betax and Fig. XandE_gamma) — illustrative of the
model's qualitative dynamics, not a forecast of actual future outcomes.

Version 2.4, May 2024
"""

import numpy as np
import matplotlib.pyplot as plt
import math

def greenGR_scenario(E0, X0, gamma0, betax, betae, alpha1, delta, mu, a, T=1000,
                      switch_t=30, betax2=None, gamma2=None):
    """
    Deterministic aggregate model (Section 3) with an optional one-time
    switch of beta_x and/or gamma at t=switch_t.

    If betax2 is None, betax stays constant for the whole simulation.
    If gamma2 is None, gamma follows the usual gamma0 - delta*Omega_t
    recursion for the whole simulation. Setting T, switch_t, betax2=None,
    gamma2=None reproduces exactly the behaviour of
    section3_baseline_model.greenGR_newnew5bis.

    Parameters
    ----------
    E0, X0 : float
        Initial net emissions and initial relative share of countries acting.
    gamma0 : float
        Baseline degree of heterogeneity (gamma with no damages), used for t < switch_t.
    betax, betae, mu : float
        Peer-pressure, private-utility, and mean-bias parameters, used for t < switch_t.
    alpha1 : float
        Effect of participation on emissions growth.
    delta : float
        Sensitivity of gamma to damages (delta=0 recovers the no-damages baseline).
    a : int
        Damage function: 0 = DICE2016R, 1 = Weitzman (2012).
    T : int, optional
        Number of periods to simulate (default 1000).
    switch_t : int, optional
        Period at which betax2/gamma2 (if given) take effect (default 30).
    betax2 : float, optional
        Value of betax used for t >= switch_t. If None, betax is unchanged.
    gamma2 : float, optional
        Value used in place of gamma0 for t >= switch_t. If None, gamma0 is unchanged.

    Returns
    -------
    E, X, v, gamma, Omega : np.ndarray
        Time series of net emissions, participation share, gamma*pi,
        degree of heterogeneity, and damages.
    """
    alpha0=0
    E=np.zeros(T)  #stock greenhouse gasses
    X=np.zeros(T)  #relative share of countries wich take climate action (-1,1)
    v=np.zeros(T)  #value of take action
    VE=np.zeros(T) #growth rate of E
    gamma=np.zeros(T) # intensity of choice 
    Omega=np.zeros(T) #damages function
    sum_E=np.zeros(T) #stock of Emissions

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
    tempAT[0]=0
    tempLO[0]=0
    Vmat[0]=3120
    Vmup[0]=5628.8
    Vmlo[0]=36706.7
    FEX[0]=0
    VFEX[0]=0.28
    VtempAT[0]=1
    VtempLO[0]=0.0068

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
    E[0]=E0
    gamma[0]=gamma0
    sum_E[0]=1730
    
    current_gamma0 = gamma0
    for t in range(T-1):
    
        if t == switch_t and gamma2 is not None:
            current_gamma0 = gamma2
            gamma[t] = gamma2
    
        current_betax = betax if (betax2 is None or t < switch_t) else betax2
    
        mat[t]=sum_E[t]+phi11*Vmat[t]+phi21*Vmup[t]
        mup[t]=phi12*Vmat[t]+phi22*Vmup[t]+phi32*Vmlo[t]
        mlo[t]=phi23*Vmat[t]+phi33*Vmlo[t]
        FEX[t]=VFEX[t]+fex
        F[t]=eta*(math.log2(mat[t]/MATP))+FEX[t]
    
        tempAT[t]=VtempAT[t]+xi1*(F[t]-xi2*VtempAT[t]-xi3*(VtempAT[t]-VtempLO[t]))
        tempLO[t]=VtempLO[t]+xi4*(VtempAT[t]-VtempLO[t])
    
        if a==0:
            Omega[t]=1-1/(1+0.0022*(tempAT[t]**2))
        else:
            Omega[t]=1-1/(1+(0.002389*(tempAT[t])**2)+(0.000005*(tempAT[t])**(6.754)))
    
        v[t]=gamma[t]*(current_betax*X[t]+betae*E[t]-mu)
        X[t+1]=((np.exp(v[t])-1)/(1+np.exp(v[t])))
        E[t+1]=E[t]+E[t]*(alpha0-alpha1*X[t])
        VE[t+1]=(E[t+1]-E[t])*100/E[t]
        gamma[t+1]=current_gamma0-delta*Omega[t]
        sum_E[t+1]=sum_E[t]+E[t]
        Vmat[t+1]=mat[t]
        Vmup[t+1]=mup[t]
        Vmlo[t+1]=mlo[t]
        VtempAT[t+1]=tempAT[t]
        VtempLO[t+1]=tempLO[t]
        VFEX[t+1]=FEX[t]

    return E, X, v, gamma, Omega

