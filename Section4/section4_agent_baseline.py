"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Agent-based model (Section 4): 2N=194 countries on a fully-connected
network, equal weights by default (weights_array=None -> uniform).
An optional weights_array parameter allows callers (e.g. the Monte
Carlo pipeline with GDP-share weights) to override this.

Version 3.0, Jan-Feb 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from network import create_network

#The network is fully connected but with weights. 
# Define the Network - the state of the Nature at the beginning of the world
#n= number of countries

# Recall the data for the weights. This are the average weight of each country, calculated as GDP_i/GDP_w, for
# the period 1989-2019. 

def GRT_AMBnetEJ(T, betan, betae, betad, N, a, mu, start_idx=0, weights_array=None):
    
    ws_graph, A, B, C=create_network(N)
    if weights_array is None:
        weights_array = np.ones(N) / N
    
    alpha0=0
    alpha1=0.02
    #T= periods
    #betan= weight of the network effect
    #betae= weight of the global effect
    #alpha0= autonomous growth rate of emission
    #alpha1= relative effect of r on E
    #N=number of countries
    i=np.ones((N,T))
    E=np.zeros(T)  #flow of net zero greenhouse gasses emissions
    X=np.zeros(T)  #relative share of countries wich take climate action (-1,1)
    sum_E=np.zeros(T) #stock of Emissions
    VE=np.zeros(T) #growth rate of E
    nC=np.zeros(T) #where I store for each t, the number of countries that take action
    
    pi=np.zeros((N, T)) # matrix when we store the v_i-values for each t in T
    y=np.zeros((N,T)) #is the matrix where, for each period, we report if countries take action (C) or abstain (A)
    Z=np.zeros((N,T)) #is the matrix for the first element of v 
    ZZ=np.zeros((N,T)) # we dived the values of Z by the degree of each variable

#initial condition: the first country acts, all others abstain
# Let's create a vector, y, size (n,1) where 1=C and 0=A and then a Matrix Z(N,T) in which we will update the single decisions
   
    y[start_idx, 0] = 1 #the starting point is the country given by start_idx (default: first country)
    Z[:,0]=A@(weights_array*y[:,0])  # in each vector there is the (weighted) number of your neighbors that take action
    
    #this matrix..should be updated...for every t in T
    ZZ[:,0]=Z[:,0]/C


    VE[0]=0
    
    nC[0]=np.sum(weights_array*y[:,0]) #calculate the total number of countries that act 
    X[0]=2*nC[0]-1
    E0=36.8
    E[0]=E0 #(alpha0-alpha1*r[0])

    #for achieving Omega we need
    D=np.zeros(T) #damages function
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
       
    #F[0]=F0
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
    sum_E[0]=1730
    error_mean = mu
    scale=1.0
    error_term = np.random.logistic(loc=error_mean, scale=scale, size=(N,T))

  
    #Loop between 0 and T for the entire model
    for t in range(T-1):
    # Damages not yet wired into the agent-based dynamics (D stays at
    # zero, betad/a have no effect) — placeholder for future work.          
    #In this new version, we do not take into account the damages' expectations 
    # but the actual value of damages. 
        mat[t]=sum_E[t]+phi11*Vmat[t]+phi21*Vmup[t]
        mup[t]=phi12*Vmat[t]+phi22*Vmup[t]+phi32*Vmlo[t]
        mlo[t]=phi23*Vmat[t]+phi33*Vmlo[t] 
        FEX[t]=VFEX[t]+fex
        F[t]=eta*(math.log2(mat[t]/MATP))+FEX[t]  # FEX exogenous ?
        
        tempAT[t]=VtempAT[t]+xi1*(F[t]-xi2*VtempAT[t]-xi3*(VtempAT[t]-VtempLO[t])) #temp atmosphere 
        tempLO[t]=VtempLO[t]+xi4*(VtempAT[t]-VtempLO[t])  #temp lower ocean

        # #When errors change over time
        # if a==0:
        #     D[t]=1-1/(1+0.0022*(tempAT[t]**2))  #omega DICE2016R
        # else: 
        #     D[t]=1-1/(1+(0.00284*(tempAT[t])**2)+(0.000005*(tempAT[t])**(6.754))) 
            
        # np.random.seed(0)  # Use the same seed for reproducibility
        # error_variance[t]=((3.1416)**2)/(3*(1-0.9*D[t])**2)
        # error_term[:,t] = np.random.logistic(loc=error_mean, scale=np.sqrt(error_variance[t]), size=N)
        
        
        
        Vmat[t+1]=mat[t]            
        Vmup[t+1]=mup[t]
        Vmlo[t+1]=mlo[t]
        VtempAT[t+1]=tempAT[t]
        VtempLO[t+1]=tempLO[t]
        VFEX[t+1]=FEX[t]
        
        #let's firstly calculate pi(:,t) for each country
        pi[:,t]=betan*(2*N*ZZ[:,t]-i[:,t])+betae*(E[t])- error_term[:,t]
   
        
        #then update counts of coutries that act and x[t]
        y[:,t+1]=np.where(pi[:,t] > 0, 1, 0) #update y if a country is C or A
        nC[t+1]=np.sum(weights_array*y[:,t+1]) #calculate the total (weighted) number of countries that act 
        
        X[t+1]=2*nC[t+1]-1  #update quota of countries that take an action, C.
        E[t+1]=E[t]+E[t]*(alpha0-alpha1*X[t]) #update emissions
        VE[t+1]=(E[t+1]-E[t])*100/E[t] #update the percentage variation of emissions
        
        sum_E[t+1]=sum_E[t]+E[t]
        #update Z, where I multiply n of neighbors that act times the inverse of degree

        Z[:,t+1]=A@(weights_array*y[:,t+1])
        ZZ[:,t+1]=Z[:,t+1]/C      
    return E, X, pi, VE, nC, y, Z, ZZ, A, C, D

np.random.seed(0)
E, X, v, VE, nC, y, Z, ZZ, A, C, D=GRT_AMBnetEJ(T=500, betan=1.1, betae=0.18, betad=0,  N=194,a=0, mu=9.77)

fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$x_t$', color = color)
plt.ylim([-1.1, 1.1])
ax1.plot(X, color = color, linewidth=1.0, linestyle='dashed', label=r'$x_t$')
ax1.tick_params(axis ='y', labelcolor = color)

ax2 = ax1.twinx()
color = 'tab:green'
ax2.set_ylabel(r'$E_t$', color = color)

ax2.plot(E, color = color,  linewidth=1.0, linestyle='dashed', label=r'E')
ax2.tick_params(axis ='y', labelcolor = color)

plt.style.use('ggplot')
plt.grid()
plt.xlim([0, 60]) 
plt.xticks([0,10,20,30,40,50,60], ['1989','1999','2009','2019','2029','2039','2049'])
plt.ylim([0,100])
fig.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2, borderaxespad=0.21)
plt.savefig('X_E_agent_baseline_singlerun.png', dpi=1200, format='png', bbox_inches='tight')
plt.show()

