# -*- coding: utf-8 -*-
"""
Created on Thu Jun  29 13:00:06 2022

@author: eim01
"""

from pre_processing_plain import *
import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

poc = 0     #point of contact
area_func = []
e_tip = 0
nu_tip = 0

fit_range = [0, 1]  #


#%%
#############################################################################
#%%
#partial unload+ mod mapping
#data import and conversion
path = 'data/Test Measurement/array-partial unload'
#path = 'data/Saphir/PM-Sa-150ms, 3um, 20nms ,3cycles'
Piezo_, MEMS_, time_, Cap_ = imp_data(path)
Piezo, MEMS, time, Cap, POC, X_val, Y_val = split_array(Piezo_, MEMS_, time_, Cap_)
#print(POC)


Results = []
indent_depth = []
fit_param = []

for i, e in enumerate(Piezo):
    
    P, M, t = data_conversion(Piezo[i], MEMS[i], time[i])
    # Piezo offset position and conversion to [nm]
    P = (P - P[0])*1000
    # detect POC
    P, M, t, poc = poc_detect(P, M, t)  
    # store each measurement as np array in a list
    Piezo[i], MEMS[i], time[i] = P, M, t
    
    #plt.plot(time[i], MEMS[i])
    
    #detect indices of load, hold and unload segment
    index, index_l, index_h, index_ul = [],[],[],[]
    index = data_splitting(P, M, t)
    index_l.append(data_splitting(P, M, t)[0])
    index_h.append(data_splitting(P, M, t)[1])
    index_ul.append(data_splitting(P, M, t)[2])
    
    #save end of segment index in list
    while index[-1] < np.argmin(M):
        ind = data_splitting(P, M, t, index_start= index[-1]+1)
        index = index + ind
        index_l.append(ind[0])
        index_h.append(ind[1])
        index_ul.append(ind[2])

    # fit curve and calculate Results
    unload_Piezo, unload_MEMS = [],[]
    popt_log, pcov_log = [], []
    S = []
    hmax = []
    
    for j in range(len(index_l)):
        up = P[index_h[j] : index_ul[j]]
        uM = M[index_h[j] : index_ul[j]]
        unload_Piezo.append(up[::-1])
        unload_MEMS.append(uM[::-1]) 
        hmax.append(up[0])
    
        
        if j==len(index_l):
            par, cov = fitting(unload_Piezo[j] , np.log(unload_MEMS[j]), [0.2, 0.95], (1.0, 1,0), fit_func=func_log)
            #uncomment for power law fit
            #par, cov = fitting(unload_Piezo[j] , np.log(unload_MEMS[j]), [0.3, 0.95], (1.0, 1.0,1), fit_func=func_exp)
        else:         
            par, cov = fitting(unload_Piezo[j] , np.log(unload_MEMS[j]), [0.3, 0.95], (1.0 ,1,0), fit_func=func_log) 
            #uncomment for power law fit
            #par, cov = fitting(unload_Piezo[j] , unload_MEMS[j], [0.3, 0.95], (1.0, 1.0,1), fit_func=func_exp)
        popt_log.append(par)
        pcov_log.append(cov)
        
       
        S.append(calc_stiff(par, up[0]))
    fit_param.append(popt_log) 
    indent_depth.append(hmax)
    Results.append(S)   #Results is a list of length(number of array points), where each entry is a list of lenght(number of load cycles)

#plot results
indent_depth = np.array(indent_depth)
Results = np.array(Results)
S_mean, S_std = [], []

for i in range(5):
    S_std.append(np.std(Results[:,i]))
    S_mean.append(np.mean(Results[:,i]))
    
#plot stiffness versus indentation depth with errorbars    
plt.errorbar(indent_depth[0], S_mean, yerr = S_std)
plt.xlabel('h [nm]')
plt.ylabel('Steifigkeit')


#     ax.plot(unload_Piezo[i], func_exp(unload_Piezo[i], par[0], par[1], par[2]), label = 'log fit' + str(i))


#%%
#from pre_processing_plain import *
path = 'data/PDMS/AR-PDMS-60-10-60-1,2um-3x3-1'   
path = 'data/PDMS/AR-PDMS-60-10-60-1,2um-4x4-1'

path = 'data/PDMS/AR-PDMS-125-10-125-2,5um-diff_LF-7x7-1'
path = 'data/PDMS/AR-PDMS-125-10-125-2,5um-diff_LF-5x5-2'
path = ['O:/5-1/5-11/Messungen/2022/06_Nico_MA/05_Datenauswertung/Python/Nanoindentation Datenauswertung/data/PDMS/AR-SA-600-1-600-4um-diff_LF-5x5-1']
path = ['data/AR-SA-200-1-200-4um-diff_LF-4x4-1']
path = ['data/AR-SA-600-1-600-4um-diff_LF-4x4-1']
path = ['data/AR-SA-600-1-600-3,5um-diff_LF-4x4-2']
path = ['data/AR-SA-335-1-335-2,5um-same_LF-5x5-2']
path = ['data/AR-SA-300-1-300-2um-diff_LF-4x4-3']
path = ['data/AR-SA-450-1-450-3um-diff_LF-4x4-1']
path = ['data/AR-SA-450-1-450-3um-diff_LF-5x5-1']
path = ['data/Saphir/AR-SA-450-1-450-3um-diff_LF-4x4-2']
path = ['data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-2']
path = ['data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-3']
#path =['O:/5-1/5-11/Messungen/2022/06_Nico_MA/05_Datenauswertung/Python/Nanoindentation Datenauswertung/data/PDMS/AR-SA-600-1-600-4um-diff_LF-5x5-1','data/AR-SA-200-1-200-4um-diff_LF-4x4-1','data/AR-SA-600-1-600-4um-diff_LF-4x4-1', 'data/AR-SA-600-1-600-3,5um-diff_LF-4x4-2', 'data/AR-SA-300-1-300-2um-diff_LF-4x4-3']
#path = ['data/AR-SA-450-1-450-3um-diff_LF-4x4-1', 'data/AR-SA-450-1-450-3um-diff_LF-5x5-1', 'data/Saphir/AR-SA-450-1-450-3um-diff_LF-4x4-2' ,'data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-2', 'data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-3']   #Messungen mit neuen MEMS
#Piezo_, MEMS_, time_, Cap_ = imp_data(path)
#Piezo, MEMS, time, Cap, POC, X_val, Y_val = split_array(Piezo_, MEMS_, time_, Cap_)
# Piezo_1, MEMS_1, time_1, Cap_1 = imp_data(path1)
# Piezo1, MEMS1, time1, Cap1, POC1, X_val1, Y_val1 = split_array(Piezo_1, MEMS_1, time_1, Cap_1)
# Piezo.extend(Piezo1)
# MEMS.extend(MEMS1)
# time.extend(time1)
# Cap.extend(Cap1)


for a, j in enumerate(path):
    Piezo_, MEMS_, time_, Cap_ = imp_data(j)
    Piezo, MEMS, time, Cap, POC, X_val, Y_val = split_array(Piezo_, MEMS_, time_, Cap_)
    S_load = []
    S_uload = []
    hmax = []
    P_in = []
    P_off = []
    #print(POC)
    n = 0
    x = 1
    if a==3:
        s=3
    else:
        s=0
    for i in range(len(Piezo)-s):
               
        P, M, t, C = data_conversion(Piezo[i*x+n], MEMS[i*x+n], time[i*x+n], Cap[i*x+n])
        # Piezo offset position and conversion to [nm]
        # detect POC
        P_, M_, t_, C_, poc = poc_detect(P, M, t, C)  
        P = (P )*1000
    
        # store each measurement as np array in a list
        Piezo[i], MEMS[i], time[i] = P, M, t
                
        Force = M * 3
        Depth = (P - M)
        Depth = Depth - Depth[np.argmin(Force[0:np.argmax(Force)])]
        # plt.subplot(2,1,2)
        # plt.plot(M, Force, label = str(i))
        # plt.grid(b=True)
        # plt.xlabel('Eindringtiefe[nm]')
        # plt.ylabel('Kraft[nN]') 
        # plt.legend()
        
        index_l, index_h, index_ul = data_splitting(P, M, t)
        
        unload_Depth = Depth[index_h : index_ul+1]
        unload_Force = Force[index_h : index_ul+1] 
    
        reversed_Depth = unload_Depth[::-1] #[nm]
        reversed_Force = unload_Force[::-1] #[nN]
        
        unload_Piezo = P[index_h : index_ul+1]
        unload_MEMS = M[index_h : index_ul+1] 
    
        reversed_Piezo = unload_Piezo[::-1] #[nm]
        reversed_MEMS = unload_MEMS[::-1] #[nN]
        
        #plot segment boundaries
        index = [index_l, index_h, index_ul]
        #calculations
        #popt_exp, pcov_exp = fitting(reversed_Piezo, reversed_MEMS, fit_range, fit_func=func_exp)
        popt_load, pcov_load = curve_fit(func_lin, P[poc:index_l], M[poc:index_l], (0.0,0.95))
        popt_uload, pcov_uload = curve_fit(func_lin, reversed_Piezo, reversed_MEMS, (0.0,0.95))
        P = P + popt_uload[1]
        S_uload.append(popt_uload[0])
        S_load.append(popt_load[0])
        hmax.append(np.max(P))
        P_in.append(Force[poc])
        P_off.append(Force[index_ul])

        #S = calc_stiff(popt_exp, reversed_Piezo[-1])    #[nN/nm]
        # h_c = calc_hc(reversed_Depth[-1], reversed_Force[-1], S, eps=0.774) #[nm]
        # E_Op, E_reduced_Op = calc_emod(S, area_sphere(h_c))
        # print(E_Op)
    
        # E_r_jkr, E_jkr = calc_JKRp(Depth, Force, R= 7500)
        # print (E_r_jkr, E_jkr)
        

        plt.subplot(3,1,1)
        plt.plot(P, M)
        plt.plot(np.take(P, index), np.take(M, index), ls = '', marker = "o")
        plt.xlabel('Piezo[nm]')
        plt.ylabel('MEMS[nm]')
        plt.grid(b=True)
        plt.title(path)  
        
    plt.subplot(3,1,2)
    plt.plot(hmax, S_load, marker = '.',label='from load segment'+str(a))
    plt.plot(hmax, S_uload, marker = '.', label='Messung'+ str(a))
    plt.grid(b=True)
    plt.xlabel('Z-max [nm]')
    plt.ylabel('relative Steigfigkeit') 
    plt.legend()
    
    plt.subplot(3,1,3)
    plt.plot(hmax, P_off, marker = '.', label ='Messung'+str(a))
    #plt.plot(hmax, P_in, marker = '.', label ='snap-in')
    plt.grid(b=True)
    plt.xlabel('Z-max [nm]')
    plt.ylabel('pull-off [nN]') 
    plt.legend()
    







