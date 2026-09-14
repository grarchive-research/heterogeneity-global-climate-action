"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Deterministic aggregate baseline model for "Heterogeneity and Global
Climate Action" (Section 3 — no agents, no network: direct iteration
of the closed-form logit map). Endogenous gamma (degree of heterogeneity)
depends on climate damages.

Version 2.4, May 2024
"""

import numpy as np 
import matplotlib.pyplot as plt
import math

 
def greenGR_newnew5bis(E0, X0, gamma0, betax, betae, alpha1, delta, mu, a, T=1000, alpha0=0.008):    
    """
    Deterministic aggregate baseline model (Section 3 of the paper).
    No agents, no network: iterates the closed-form logit map directly.

    Equations:
        pi_t = betax*X_t + betae*E_t - mu                (eq. 3/8)
        X_{t+1} = (exp(gamma_t*pi_t) - 1) / (1 + exp(gamma_t*pi_t))   (eq. 11/14)
        E_{t+1} = E_t * (1 + alpha0 - alpha1*X_t)         (eq. 9), alpha0 fixed at 0
        gamma_{t+1} = gamma0 - delta*Omega_t              (eq. 16, endogenous heterogeneity)
        Omega_t = damage function (DICE if a==0, Weitzman otherwise), driven by
                  the carbon-cycle block (temperature, radiative forcing).

    Parameters
    ----------
    E0, X0 : float
        Initial net emissions and initial relative share of countries acting.
    gamma0 : float
        Baseline degree of heterogeneity (gamma with no damages).
    betax, betae, mu : float
        Peer-pressure, private-utility, and mean-bias parameters.
    alpha1 : float
        Effect of participation on emissions growth.
    delta : float
        Sensitivity of gamma to damages (delta=0 recovers the no-damages baseline).
    a : int
        Damage function: 0 = DICE2016R, 1 = Weitzman (2012).
    T : int, optional
        Number of periods to simulate (default 1000).

    Returns
    -------
    E, X, v, gamma, Omega : np.ndarray
        Time series of net emissions, participation share, gamma*pi, 
        degree of heterogeneity, and damages.
    """
    
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
    for t in range(T-1):  
        
        
        mat[t]=sum_E[t]+phi11*Vmat[t]+phi21*Vmup[t]
        mup[t]=phi12*Vmat[t]+phi22*Vmup[t]+phi32*Vmlo[t]
        mlo[t]=phi23*Vmat[t]+phi33*Vmlo[t] 
        FEX[t]=VFEX[t]+fex
        F[t]=eta*(math.log2(mat[t]/MATP))+FEX[t] 
        
        tempAT[t]=VtempAT[t]+xi1*(F[t]-xi2*VtempAT[t]-xi3*(VtempAT[t]-VtempLO[t])) #temp atmosphere 
        tempLO[t]=VtempLO[t]+xi4*(VtempAT[t]-VtempLO[t])  #temp lower ocean
        

        if a==0:
            Omega[t]=1-1/(1+0.0022*(tempAT[t]**2))  #omega DICE2016R
        else: 
            Omega[t]=1-1/(1+(0.002389*(tempAT[t])**2)+(0.000005*(tempAT[t])**(6.754))) 

    
        
        v[t]=gamma[t]*(betax*X[t]+betae*E[t]-mu) 
        X[t+1]=((np.exp(v[t])-1)/(1+np.exp(v[t])))
        E[t+1]=E[t]+E[t]*(alpha0-alpha1*X[t])
        VE[t+1]=(E[t+1]-E[t])*100/E[t]
        gamma[t+1]=gamma0-delta*Omega[t]
        sum_E[t+1]=sum_E[t]+E[t]
        Vmat[t+1]=mat[t]            
        Vmup[t+1]=mup[t]
        Vmlo[t+1]=mlo[t]
        VtempAT[t+1]=tempAT[t]
        VtempLO[t+1]=tempLO[t]
        VFEX[t+1]=FEX[t]
  
    return E, X, v, gamma, Omega

if __name__ == "__main__":

    E, X, v, gamma, Omega= greenGR_newnew5bis(E0=36.8, X0=-.97, gamma0=1, betax=1.1, betae=0.18, alpha1=0.014, delta=0, mu=9.77, a=0)
    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.set_xlabel(r'$t$')
    ax1.set_ylabel(r'$x_t$', color = color)
    plt.ylim([-1.1, 1.1])
    plt.xlim([0, 100])
    ax1.plot(X, color = color, linewidth=1.0, linestyle='dashed', label=r'$x_t$')
    ax1.tick_params(axis ='y', labelcolor = color)

    ax2 = ax1.twinx()
    color = 'tab:green'
    ax2.set_ylabel(r'$E_t$', color = color)

    ax2.plot(E, color = color,  linewidth=1.0, linestyle='dashed', label=r'$E_t$')
    ax2.tick_params(axis ='y', labelcolor = color)

    plt.style.use('ggplot')
    plt.grid()
    plt.xlim([0, 60])
    plt.xticks([0,10,20,30,40,50,60], ['1989','1999','2009','2019','2029','2039','2049'])
    plt.ylim([0,100])

    fig.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2, borderaxespad=0.21)
    plt.savefig('X_E_estimatedparameters_alpha0.png', dpi=1200, format='png', bbox_inches='tight') 
    plt.show()