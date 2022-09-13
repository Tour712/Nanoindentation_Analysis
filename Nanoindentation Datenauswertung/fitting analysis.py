# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 10:41:47 2022

@author: eim01
"""
from functions import *
import numpy as np
#evaluation of fitting algorithm

iterations = 200
start_param = (0.05,1.25,0)
h_max  = 200    #[nm]
n_points = [500]
noise_std_P = np.linspace(0,4,10)
alpha = 0.05    #in [mN/nm^m]
m = 1.25
h_f = 0
h = np.linspace(30, h_max, 500)#[nm]
P_opt = alpha*(h - h_f)**m 
#P_t = (0.1)*(h + 10)**(1.1)
S = calc_stiff([alpha, m , h_f], h_max)

plt.plot (h, P_opt)
#plt.plot (h, P_n)
plt.title('ideale Daten')
plt.xlabel('h in [nm]')
plt.ylabel('P in [mN]')   


#popt_exp_g, cov_exp_g = [], []
popt_exp, cov_exp = [], []
popt_std_alpha, popt_std_m, popt_std_hf = [], [], []
popt_mean_alpha, popt_mean_m, popt_mean_hf = [], [], []
S_std = []
S_mean = []
#%%

for k in n_points:
    popt_std_alpha, popt_std_m, popt_std_hf = [], [], []
    S_std = []
    h = np.linspace(0, h_max, k)#[nm]
    P_opt = alpha*(h - h_f)**m
    for n,i in enumerate(noise_std_P):
        h = np.linspace(0, h_max, k)#[nm]
        popt = np.empty([iterations,3])
        S_t = np.empty([iterations,1])
        
        for j in range(iterations):
            noise_P = np.random.normal(0, 1, k)
            P = P_opt + (P_opt[-1]*(i/100)) *noise_P    #[mN] prozentual
            #P = P_opt + noise_P*0.001    #[mN] absolut
            p, c = fitting(h, P, fit_range, start_param, fit_func=func_exp)
            popt[j] = p
            S_t[j] = calc_stiff(popt[j], h_max)
        
        plt.subplot(int(len(noise_std_P)/2), int(len(noise_std_P)/2),n+1)
        plt.plot(h,P)
        
        popt_std_alpha.append((np.std(popt[:,0])/alpha)*100)
        popt_std_m.append((np.std(popt[:,1])/m)*100)
        popt_std_hf.append(np.std(popt[:,2]))
        popt_mean_alpha.append(np.mean(popt[:,0]))
        popt_mean_m.append(np.mean(popt[:,1]))
        popt_mean_hf.append(np.mean(popt[:,2]))
        S_std.append((np.std(S_t)/S)*100) 
        #S_std.append(np.std(S_t)) 
        S_mean.append(np.mean(S_t))
    
    plt.figure()
    plt.subplot(2,2,1)
    plt.plot(noise_std_P, popt_std_alpha, label=str(k)+' datapoints')
    plt.title('standard deviation of alpha as a function of noise level (alpha=0,05)')
    plt.xlabel('noise level in [%]')
    plt.ylabel('standard deviation in [%]')   
    plt.legend()
    
    plt.subplot(2,2,2)
    plt.plot(noise_std_P, popt_std_m, label=str(k)+' datapoints')
    plt.title('standard deviation of m as a function of noise level (m=1.25)')
    plt.xlabel('noise level in [%]')
    plt.ylabel('standard deviation in [%]')     
    plt.legend()
    
    plt.subplot(2,2,3)
    plt.plot(noise_std_P, popt_std_hf, label=str(k)+' datapoints')
    plt.title('standard deviation of hf as a function of noise level (hf=0)')
    plt.xlabel('noise level [nm]')
    plt.ylabel('standard deviation [nm]') 
    plt.legend()  
    
    plt.subplot(2,2,4)
    plt.plot(noise_std_P, S_std, label=str(k)+' datapoints')
    plt.title('standard deviation of S as a function of noise level (S=0.235)')
    plt.xlabel('noise level [%]')
    plt.ylabel('standard deviation [%]') 
    plt.legend() 
    
     
    #popt_exp_g.append()
    #a,m_t,h_t = 0,0,0
    #calculate mean values (doesnt make sense)
    # for i in popt_exp:
    #     a += i[0]
    #     m_t += i[1]
    #     h_t += i[2]    
    # a = a/len(popt_exp)    
    # m_t = m_t/len(popt_exp)
    # h_f = h_f/len(popt_exp)
    #popt_exp_g.append([a,m_t,h_f])    