import numpy as np
import math


def g(corr,bin_edges,nr_steps,L,act_N):
    delta_R=bin_edges[-1]/(len(bin_edges)-1)
    R=bin_edges[0:len(bin_edges)-1] + delta_R/2
    av_corr= np.divide(np.mean(corr,1),R*R)/(4*math.pi*delta_R*act_N*(act_N-1))*L**3
    sigma_corr=np.zeros(len(bin_edges)-1)
    
    #Bootstrapping error estimation
    sample=np.zeros(500) #Resampling 500 times
    for j in range(0,len(bin_edges)-1):#Sampling for each bin
        for i in range(0,500): 
            sample[i]=np.mean(np.random.choice(corr[j,:], nr_steps)) #For  each bin, we take the mean of nr_steps randomly chosen values
        sigma_corr[j]=np.sqrt(np.var(sample))
    sigma_av_corr= np.divide(sigma_corr,R*R)/(4*math.pi*delta_R*act_N*(act_N-1))*L**3
    return av_corr, R, sigma_corr

def p(press, nr_steps):
    
    mean_press=np.mean(press)
    
    #Bootstrapping error estimation
    sample=np.zeros(500)
    for i in range(0,500): #Resampling 500 times
            sample[i]=np.mean(np.random.choice(press, nr_steps)) #For  each bin, we take the mean of nr_steps randomly chosen values
    sigma_press=np.sqrt(np.var(sample))
    return mean_press,  sigma_press
