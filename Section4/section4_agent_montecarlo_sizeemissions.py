"""
@author: Giorgio Ricchiuti
@website: www.grarchive.net

Description
-----------
Monte Carlo pipeline for the agent-based model (Section 4), GDP-share
weights and country-level emissions: 2N=194 countries on a fully-connected
network. Each country's own emissions Ei grow or shrink by a fixed rate
depending on whether it acts, and global net emissions E are the sum
across countries. Nested double Monte Carlo — idiosyncratic error term
(inner loop) x 194 starting countries (outer loop). Produces the
published Fig. Comparison_evolution_discrete_size_emissions.

Version 3.0, Jan-Feb 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from dataset import dataset
from network import create_network

weights_array, GHG_array = dataset()

def MonteCarlo3(N, T, MC_runs):
    
    E_results = np.zeros((N, MC_runs, T))  # Stores emissions for each country and MC run
    X_results = np.zeros((N, MC_runs, T))  # Stores action shares
    VE_results = np.zeros((N, MC_runs, T))  # Stores emission growth rates
    nC_results = np.zeros((N, MC_runs, T))  # Stores number of acting countries
    
    for k in range(N):  # Iterate over each country as the starting point
        
        for mc in range(MC_runs):  # Monte Carlo loop for error term variations
            np.random.seed(mc)  #set the seed
            
            def GRT_AMBnetEJ_EiMC(T, betan, betae, betad, N, a, mu):
                
                ws_graph, A, B, C=create_network(N)
                
                #alpha0=0 #autonomous growth rate of emission
                #alpha1=0.02 #relative effect of x on E
                #T= periods
                #betan= weight of the network effect
                #betae= weight of the global effect
                #N=number of countries
                i=np.ones((N,T))
                E=np.zeros(T)  #flow of net zero greenhouse gasses emissions
                Ei=np.zeros((N,T)) #flow of net zero GHG at country level
                X=np.zeros(T)  #relative share of countries wich take climate action (-1,1)
                sum_E=np.zeros(T) #stock of Emissions
                VE=np.zeros(T) #growth rate of E
                nC=np.zeros(T) #where I store for each t, the number of countries that take action
                
                #to get the the matrix v
                #Ei=np.zeros((N,T)) # flows of net zero emissions at country level
                pi=np.zeros((N, T)) # matrix when we store the p_i-values for each t in T
                y=np.zeros((N,T)) #is the matrix where, for each period, we report if countries take action (C) or abstain (A)
                Z=np.zeros((N,T)) #is the matrix for the first element of v 
                ZZ=np.zeros((N,T)) # we dived the values of Z by the degree of each variable
                
                g=np.zeros((N,T)) #growth rate at country level due to action or not

                #initial conditions
                #Initialization: Randomly we decide if countries take action (C) or abstain (A)
                # Let's create a vector, y, size (n,1) where 1=C and 0=A and then a Matrix Z(N,T) in which we will update the single decisions
               
               
                y[k, 0] = 1
                #y[:,0]=np.random.choice([0, 1], size=N) #this is a vector variable 
                #weights_array=np.ones(N) / N
                Z[:,0]=A@(weights_array*y[:,0])  # in each vector there is the (weighted) number of your neighbors that take action
                
                #this matrix..should be updated...for every t in T
                ZZ[:,0]=Z[:,0]/C

                # The initial conditions here will be the same, as for the other paper. 
                VE[0]=0
                
                nC[0]=np.sum(weights_array*y[:,0]) #calculate the total number of countries that act 
                X[0]=2*nC[0]-1
                Ei[:,0]=GHG_array
                
                E0=np.sum(Ei[:, 0])  #It should be 36.5 very close to 36.8
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
                sum_E[0]=1730
                # Idiosyncratic logistic error term
                
                error_mean = mu
                #error_variance=3.29
                #xerror_variance=1.047
                scale=1.0
                error_term=np.zeros((N,T))
                #error_variance=np.zeros(T)
                error_term = np.random.logistic(loc=error_mean, scale=scale, size=(N,T))
                
                phi=0.02
                varphi=0.02
                
              
                #Loop between 0 and T for the entire model
                for t in range(T-1):
                      
                #In this version, we do not take into account the damages' expectations 
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
                    
                    #let's firstly calculate v(:,t) for each country
                    #pi[:,t]=betan*(2*N*ZZ[:,t]-i[:,t])+betae*(Ei[:,t])- error_term[:,t] #with have country size and emissions
                    
                    pi[:,t]=betan*(2*N*ZZ[:,t]-i[:,t])+betae*(E[t])- error_term[:,t] #with have country size and emissions
                    g[:,t] = np.where(pi[:, t] > 0, -varphi, phi)  # Compute g^i_t based on action
                    Ei[:, t+1] = Ei[:, t] * (1 + g[:,t])  # Apply the update equation
                    E[t+1] = np.sum(Ei[:, t+1])  # Compute total emissions at time t
                    #then update counts of coutries that act and x[t]
                    y[:,t+1]=np.where(pi[:,t] > 0, 1, 0) #update y if a country is C or A
                    nC[t+1]=np.sum(weights_array*y[:,t+1]) #calculate the total (weighted) number of countries that act 
                    
                    X[t+1]=2*nC[t+1]-1  #update quota of countries that take an action, C.
                    #E[t+1]=E[t]+E[t]*(alpha0-alpha1*X[t]) #update emissions
                            
                    VE[t+1]=(E[t+1]-E[t])*100/E[t] #update the percentage variation of emissions
                    sum_E[t+1]=sum_E[t]+E[t]
               
                    Z[:,t+1]=A@(weights_array*y[:,t+1])
                    ZZ[:,t+1]=Z[:,t+1]/C      
                    #Let's close the loop. Return different variables. 

                return E, Ei, X, pi, VE, nC, y, Z, ZZ, A, C, D
        
        
            print(f"Running model with initial country {k}")
        
            
            # Run the model with initial condition at index k
            E, Ei, X, v, VE, nC, y, Z, ZZ, A, C, D=GRT_AMBnetEJ_EiMC(T=500, betan=1.1, betae=0.18, betad=0,  N=194,a=0, mu=9.77)
        

   
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

mean_E, mean_VE, mean_X, mean_nC, std_E, std_X, std_VE, std_nC= MonteCarlo3(194, 500, MC_runs=50)

# Styling
plt.style.use('ggplot')
plt.grid()


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

plt.xlim([0, 60]) 
plt.xticks([0,10,20,30,40,50,60], ['1989','1999','2009','2019','2029','2039','2049'])
plt.ylim([20, 80])
fig.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2, borderaxespad=0.21)
plt.savefig('MonteCarlo_withcountriesweights_emissions.png', dpi=1200, format='png', bbox_inches='tight')
plt.show()
