"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Monte Carlo pipeline for the agent-based model (Section 4) comparing
three damage-function assumptions — no damages (baseline), DICE2016R,
and Weitzman (2012) — with a one-time peer-pressure switch at t=30:
2N=194 countries on a fully-connected network, equal weights. Nested
double Monte Carlo — 100 draws of the idiosyncratic error term (inner
loop) x 194 starting countries (outer loop). Produces the published
Fig. Comparison_IAM_evolution_discrete.

Version 3.0, Jan-Feb 2025
"""


import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
import math
from network import create_network

# Define the Network - the state of the Nature at the beginning of the world
#n= number of countries



#Let's start the model. 


def MonteCarlo_200(N, T, MC_runs, a, b):
    
    E_results = np.zeros((N, MC_runs, T))  # Stores emissions for each country and MC run
    X_results = np.zeros((N, MC_runs, T))  # Stores action shares
    VE_results = np.zeros((N, MC_runs, T))  # Stores emission growth rates
    nC_results = np.zeros((N, MC_runs, T))  # Stores number of acting countries
    
    for k in range(N):  # Iterate over each country as the starting point
        
        for mc in range(MC_runs):  # Monte Carlo loop for error term variations
            np.random.seed(mc)  #set the seed


            def GRT_AMBnetIAM(T, betan, betan2, betae, betad, N,a, mu, b):
                # Parameters:
                # - N (int): Number of nodes in the network.
                # - k (int): Number of nearest neighbors to connect each node to.
                # - p (float): Probability of rewiring each edge.
                #     Returns:
                # - ws_graph (networkx.Graph): The generated small-world network.
                # - A (numpy.ndarray): Adjacency matrix of the network.
                # - B (numpy.ndarray): Degree array of nodes.
                # - C (numpy.ndarray): Degrees of nodes.
                # - A_with_diagonal (networkx.Graph): Graph with self-loops added.
                #Maybe we should insert the network's generation within the model. Let'see
                #ws_graph, A, C,D=Network(n,2,0)
                # a=0 DICE model, a=1 Weitzman
                ws_graph, A, B, C=create_network(N)
                
                alpha0=0
                alpha1=0.02
                #T= periods
                #betan= weight of the network effect
                #betaG= weight of the global effect
                #betad= weight of the climate damages
                #alpha0= autonomous growth rate of emission
                #alpha1= relative effect of r on E
                #N=number of countries
                #M= maximun periods of forward looking, calculate the Expectation of Omega after 10,20,30 periods. 
                i=np.ones((N,T))
                E=np.zeros(T)  #stock greenhouse gasses
                X=np.zeros(T)  #relative share of countries wich take climate action (-1,1)
                sum_E=np.zeros(T) #stock of Emissions
                v=np.zeros(T)  #value of take action
                VE=np.zeros(T) #growth rate of E
                nC=np.zeros(T) #where I store for each t, the number of countries that take action
                
                #to get the the matrix v
                v=np.zeros((N, T)) # matrix when we store the v_i-values for each t in T
                y=np.zeros((N,T)) #is the matrix where, for each period, we report if countries take action (C) or abstain (A)
                Z=np.zeros((N,T)) #is the matrix for the first element of v 
                ZZ=np.zeros((N,T)) # we dived the values of Z by the degree of each variable
            
                
            #initial conditions
            #Initialization: Randomly we decide if countries take action (C) or abstain (A)
            # Let's create a vector, y, size (n,1) where 1=C and 0=A and then a Matrix Z(N,T) in which we will update the single decisions
               
                #np.random.seed(0)  #set the seed
                y[k, 0] = 1
                #y[:,0]=np.random.choice([0, 1], size=N) #this is a vector variable 
                
                #I=np.identity(N) #we need this for the loop. 
                
                Z[:,0]=A@y[:,0]  # in each vector there is the number of your neighbors that take action
                
                #this matrix..should be updated...for every t in T
                ZZ[:,0]=Z[:,0]/C
            
                # The initial conditions here will be the same, as for the other paper. 
                VE[0]=0
                
                nC[0]=np.sum(y[:,0]) #calculate the total number of countries that act 
                X[0]=(2*nC[0]-N)/N 
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
                # Add a logistically distributed error term
            
                # if b==0:
                #     np.random.seed(0)  # Use the same seed for reproducibility
                #     error_mean = mu
                error_term=np.zeros((N,T))
                #     error_variance=np.zeros(T)
                
            
                error_mean = mu
                delta=0.99  
            
                
                #Loop between 0 and T for the entire model
                for t in range(30):
                      
                #In this new version, we do not take into account the damages' expectations 
                # but the actual value of damages. 
                    mat[t]=sum_E[t]+phi11*Vmat[t]+phi21*Vmup[t]
                    mup[t]=phi12*Vmat[t]+phi22*Vmup[t]+phi32*Vmlo[t]
                    mlo[t]=phi23*Vmat[t]+phi33*Vmlo[t] 
                    FEX[t]=VFEX[t]+fex
                    F[t]=eta*(math.log2(mat[t]/MATP))+FEX[t]  # FEX exogenous ?
                    
                    tempAT[t]=VtempAT[t]+xi1*(F[t]-xi2*VtempAT[t]-xi3*(VtempAT[t]-VtempLO[t])) #temp atmosphere 
                    tempLO[t]=VtempLO[t]+xi4*(VtempAT[t]-VtempLO[t])  #temp lower ocean
            
                    if b==0:
                    
                        if a==0:
                            D[t]=1-1/(1+0.0022*(tempAT[t]**2))  #omega DICE2016R
                        else: 
                            D[t]=1-1/(1+(0.002389*(tempAT[t])**2)+(0.000005*(tempAT[t])**(6.754))) 
                            
                        np.random.seed(mc)  # Use the same seed for reproducibility
                        scale2=1/(1-delta*D[t])
                        #error_term = np.random.logistic(loc=error_mean, scale=scale, size=(N,T))
                        error_term = np.random.logistic(loc=error_mean, scale=scale2, size=(N,T))
                    
                    else: 
                        np.random.seed(mc)  # Use the same seed for reproducibility
                        scale=1
                        error_term = np.random.logistic(loc=error_mean, scale=scale, size=(N,T))
                    
                    
                    Vmat[t+1]=mat[t]            
                    Vmup[t+1]=mup[t]
                    Vmlo[t+1]=mlo[t]
                    VtempAT[t+1]=tempAT[t]
                    VtempLO[t+1]=tempLO[t]
                    VFEX[t+1]=FEX[t]
                    
                    #let's firstly calculate v(:,t) for each country
                    v[:,t]=betan*(2*ZZ[:,t]-i[:,t])+betae*(E[t])- error_term[:,t]
               
                    
                    #then update counts of coutries that act and x[t]
                    y[:,t+1]=np.where(v[:,t] > 0, 1, 0) #update y if a country is C or A
                    nC[t+1]=np.sum(y[:,t+1]) #calculate the total number of countries that act 
                    
                    X[t+1]=(2*nC[t+1]-N)/N  #update quota of countries that take an action, C.
                    E[t+1]=E[t]+E[t]*(alpha0-alpha1*X[t]) #update emissions
                    VE[t+1]=(E[t+1]-E[t])*100/E[t] #update the percentage variation of emissions
                    
                    sum_E[t+1]=sum_E[t]+E[t]
                    #update Z, where I multiply n of neighbors that act times the inverse of degree
                    #Y[t+1]=I*y[:,t+1]
                    #Z[:,t+1]=np.dot(Y[t+1],D).flatten()            
                    Z[:,t+1]=A@y[:,t+1]
                    ZZ[:,t+1]=Z[:,t+1]/C      
                    #Let's close the loop. Return different variables. 
            
                
                
                for t in range(30,T-1):
                
                    #In this new version, we do not take into account the damages' expectations 
                    # but the actual value of damages. 
                        mat[t]=sum_E[t]+phi11*Vmat[t]+phi21*Vmup[t]
                        mup[t]=phi12*Vmat[t]+phi22*Vmup[t]+phi32*Vmlo[t]
                        mlo[t]=phi23*Vmat[t]+phi33*Vmlo[t] 
                        FEX[t]=VFEX[t]+fex
                        F[t]=eta*(math.log2(mat[t]/MATP))+FEX[t]  # FEX exogenous ?
                        
                        tempAT[t]=VtempAT[t]+xi1*(F[t]-xi2*VtempAT[t]-xi3*(VtempAT[t]-VtempLO[t])) #temp atmosphere 
                        tempLO[t]=VtempLO[t]+xi4*(VtempAT[t]-VtempLO[t])  #temp lower ocean
            
                        # When errors change over time
                        
                        if b==0:
                                        
                            if a==0:
                                D[t]=1-1/(1+0.0022*(tempAT[t]**2))  #omega DICE2016R
                            else: 
                                D[t]=1-1/(1+(0.002389*(tempAT[t])**2)+(0.000005*(tempAT[t])**(6.754)))                                
                            #np.random.seed(0)  # Use the same seed for reproducibility
                            scale2=1/(1-delta*D[t])
                            error_term= np.random.logistic(loc=error_mean, scale=scale2,  size=(N,T))                
                        else: 
                            #np.random.seed(0)  # Use the same seed for reproducibility
                            scale=2
                            error_term = np.random.logistic(loc=error_mean, scale=scale, size=(N,T))
                  
                        
                        Vmat[t+1]=mat[t]            
                        Vmup[t+1]=mup[t]
                        Vmlo[t+1]=mlo[t]
                        VtempAT[t+1]=tempAT[t]
                        VtempLO[t+1]=tempLO[t]
                        VFEX[t+1]=FEX[t]
                        
                        #let's firstly calculate v(:,t) for each country
                        #v[:,t]=betan2*(2*ZZ[:,t]-i[:,t])+betae*(E[t])- error_term
            
                        v[:,t]=betan2*(2*ZZ[:,t]-i[:,t])+betae*(E[t])- error_term[:,t]
                   
                        
                        #then update counts of coutries that act and x[t]
                        y[:,t+1]=np.where(v[:,t] > 0, 1, 0) #update y if a country is C or A
                        nC[t+1]=np.sum(y[:,t+1]) #calculate the total number of countries that act 
                        
                        X[t+1]=(2*nC[t+1]-N)/N  #update quota of countries that take an action, C.
                        E[t+1]=E[t]+E[t]*(alpha0-alpha1*X[t]) #update emissions
                        VE[t+1]=(E[t+1]-E[t])*100/E[t] #update the percentage variation of emissions
                        
                        sum_E[t+1]=sum_E[t]+E[t]
                        #update Z, where I multiply n of neighbors that act times the inverse of degree
                     
                        Z[:,t+1]=A@y[:,t+1]
                        ZZ[:,t+1]=Z[:,t+1]/C      
                        #Let's close the loop. Return different variables. 
                
                return E, X, v, VE, nC, y, Z, ZZ, A
            
            print(f"Running model with initial country {k}")


            E, X, v, VE, nC, y, Z, ZZ, A=GRT_AMBnetIAM(T=100, betan=1.1, betan2=1.1, betae=0.18, betad=0,  N=194, a=a, mu=9.77, b=b)

   
            # Store results for this MC iteration
            #print(f"E shape: {E.shape}, Expected shape: {E_results[k, mc, :].shape}")
            E_results[k, mc, :] = E
            X_results[k, mc, :] = X
            VE_results[k, mc, :] = VE
            nC_results[k, mc, :] = nC
            
        # Compute the mean over all MC runs
        mean_E = np.mean(E_results, axis=(0,1))  # Averaging over both N and MC runs
        mean_X = np.mean(X_results, axis=(0,1))
        mean_VE = np.mean(VE_results, axis=(0,1))
        mean_nC = np.mean(nC_results, axis=(0,1))
        
        # Compute standard deviation over all MC runs
        std_E = np.std(E_results, axis=(0,1))
        std_X = np.std(X_results, axis=(0,1))
        std_VE = np.std(VE_results, axis=(0,1))
        std_nC = np.std(nC_results, axis=(0,1))
        
    return mean_E, mean_VE, mean_X, mean_nC, std_E, std_X, std_VE, std_nC


mean_E1, _, mean_X1, _, std_E1, std_X1, _, _ = MonteCarlo_200(194, 100, MC_runs=50, a=0,b=1)
mean_E2, _, mean_X2, _, std_E2, std_X2, _, _ = MonteCarlo_200(194, 100, MC_runs=50, a=1,b=0)
mean_E3, _, mean_X3, _, std_E3, std_X3, _, _ = MonteCarlo_200(194, 100, MC_runs=50, a=0,b=0)

T = 100

# --- First Figure: X_t ---
fig1, ax1 = plt.subplots(figsize=(8, 5))  # Set figure size
ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$X_t$', color='red')
ax1.tick_params(axis='y', labelcolor='red')
# Plot mean X_t for all cases
ax1.plot(mean_X1, color='tab:red', linewidth=1.0, linestyle='dashed', label=r'Baseline')
ax1.fill_between(range(T), mean_X1 - std_X1, mean_X1 + std_X1, color='tab:red', alpha=0.1)
ax1.plot(mean_X2, color='tab:blue', linewidth=1.0, linestyle='dashed', label=r'Weitzman')
ax1.fill_between(range(T), mean_X2 - std_X2, mean_X2 + std_X2, color='tab:blue', alpha=0.1)
ax1.plot(mean_X3, color='black', linewidth=1.0, linestyle='dotted', label=r'DICE')
ax1.fill_between(range(T), mean_X3 - std_X3, mean_X3 + std_X3, color='black', alpha=0.1)  # Fixed color

ax1.set_ylim([-1.01, 1.01]) 
ax1.set_xlim([0, 60]) 
plt.xticks([0,10,20,30,40,50,60], ['1989','1999','2009','2019','2029','2039','2049'])
ax1.set_xlim([0, 60])
ax1.grid(True)
ax1.legend(loc='lower right', fontsize=8)
plt.savefig('Comparison_IAM_evolution_discrete.png', dpi=1200, format='png', bbox_inches='tight')
plt.show()
