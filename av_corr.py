#This function requires the of corr and its calculation for every time step, as well as bin_edges


import numpy as np
import math


def g(corr,bin_edges,nr_steps,L,act_N):
    delta_R=bin_edges[-1]/(len(bin_edges)-1)
    R=bin_edges[0:len(bin_edges)-1] + delta_R/2
    av_corr= np.divide(np.mean(corr,1),R*R)/(4*math.pi*delta_R*act_N*(act_N-1))*L**3
    
    #Bootstrapping error
    sigma_corr=np.zeros(len(bin_edges)-1)
    for j in range(0,len(bin_edges)-1):#Sampling for each bin
        for i in range(0,500): #Resampling 500 times
            sample=np.random.choice(corr[j,:], nr_steps) #For  each bin, we take as many samples as time steps
        sigma_corr[j]=np.sqrt(np.var(sample)/(nr_steps-1))
    
    return av_corr, R, sigma_corr
