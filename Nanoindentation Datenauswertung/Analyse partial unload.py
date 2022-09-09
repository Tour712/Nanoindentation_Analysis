# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 10:24:23 2022

@author: eim01
"""

from pre_processing_plain import *
###########################################################################
#this script is for analysis of partial-unload measurement on Saphire

path = 'data/Saphir/PM-Sa-150ms, 3um, 20nms ,9cycles-3'
#path = 'data/Saphir/PM-Sa-150ms, 3um, 7nms ,3cycles-2'


Piezo, MEMS, time, Cap = imp_data(path)
Piezo_np, MEMS_np, time_np, Cap_np = data_conversion(Piezo, MEMS, time, Cap)
Piezo_np = (Piezo_np - Piezo_np[0])*1000

Force = MEMS_np*K

P, M, t, C, poc_i = poc_detect(Piezo_np, MEMS_np, time_np, Cap_np)

index, index_l, index_h, index_ul = [],[],[],[]

index = data_splitting(Piezo_np, MEMS_np, time_np)
index_l.append(data_splitting(Piezo_np, MEMS_np, time_np)[0])
index_h.append(data_splitting(Piezo_np, MEMS_np, time_np)[1])
index_ul.append(data_splitting(Piezo_np, MEMS_np, time_np)[2]) 

# #save end of segment index in list
while index[-1] < np.argmin(MEMS_np):
    ind = []
    ind = data_splitting(Piezo_np, MEMS_np, time_np, index_start= index[-1]+1)
    index = index + ind
    index_l.append(ind[0])
    index_h.append(ind[1])
    index_ul.append(ind[2])
    
plt.subplot(3,1,1)
plt.plot(time_np ,MEMS_np, ls = '', marker = "+", markersize = 2)
plt.plot(np.take(time_np, index), np.take(MEMS_np, index), ls = '', marker = "o", label = 'segment boundarys')
plt.xlabel('Zeit [s]')
plt.ylabel('MEMS Verschiebung [nm]')
plt.legend()

plt.subplot(3,1,2)
plt.plot(Piezo_np, MEMS_np)
plt.plot(np.take(Piezo_np, index), np.take(MEMS_np, index), ls = '', marker = "o", label = 'segment boundarys')
plt.xlabel('Piezoposition [nm]')
plt.ylabel('MEMS Verschiebung [nm]')


unload_Piezo, unload_MEMS = [],[]
popt_log, pcov_log = [],[]
S = []
hmax = []
for i in range(len(index_l)):
    up = Piezo_np[index_h[i] : index_ul[i]]
    uM = MEMS_np[index_h[i] : index_ul[i]]
    unload_Piezo.append(up[::-1])
    unload_MEMS.append(uM[::-1]) 
    hmax.append(np.max(up))
    if i==len(index_l):
        par, cov = fitting(unload_Piezo[i] , unload_MEMS[i], [0.2, 0.95])
        #uncomment for power law fit
        #par, cov = fitting(unload_Piezo[i] , np.log(unload_MEMS[i]), [0.4, 0.95], (1.0, 1.0, 95), fit_func=func_exp)
    else:         
        par, cov = fitting(unload_Piezo[i] , unload_MEMS[i], [0.3, 0.95]) 
        #par, cov = fitting(unload_Piezo[i] , unload_MEMS[i], [0.3, 0.95], (1.0 ,1,0)) 
        #uncomment for power law fit
        #par, cov = fitting(unload_Piezo[i] , unload_MEMS[i], [0.4, 0.95], (1.0, 1.0 , 95), fit_func=func_exp)
    popt_log.append(par)
    pcov_log.append(cov)
    S.append(par[0])
    
plt.subplot(3,1,3)
plt.plot(hmax, S, marker = '.')

#######################################################################################################
#%%
#This script is for the analysis of an array measurement on sapphire, with partial unload load function

path = 'data/Saphir/APM-Sa-150ms, 3um, 7nms, 16 cycles-4'

Piezo_, MEMS_, time_, Cap_ = imp_data(path)
Piezo, MEMS, time, Cap, POC, X_val, Y_val = split_array(Piezo_, MEMS_, time_, Cap_)
S_load = []
S_uload = []
hmax = []
P_in = []
P_off = []


for j in range(len(Piezo)):
    
    P, M, t, C = data_conversion(Piezo[j*x+n], MEMS[j*x+n], time[j*x+n], Cap[j*x+n])
    # Piezo offset position and conversion to [nm]
    # detect POC
    P_, M_, t_, C_, poc = poc_detect(P, M, t, C)  
    P = (P-P[0] )*1000

    # store each measurement as np array in a list
    Piezo[j], MEMS[j], time[j] = P, M, t
       
    index, index_l, index_h, index_ul = [],[],[],[]
    
    index = data_splitting(Piezo_np, MEMS_np, time_np)
    index_l.append(data_splitting(P, M, t)[0])
    index_h.append(data_splitting(P, M, t)[1])
    index_ul.append(data_splitting(P, M, t)[2]) 
    
    # #save end of segment index in list
    while index[-1] < np.argmin(M):
        ind = []
        ind = data_splitting(P, M, t, index_start= index[-1]+1)
        index = index + ind
        index_l.append(ind[0])
        index_h.append(ind[1])
        index_ul.append(ind[2])
        
    # plt.subplot(3,1,1)
    # plt.plot(t ,M, ls = '', marker = "+", markersize = 2)
    # plt.plot(np.take(t, index), np.take(M, index), ls = '', marker = "o", label = 'segment boundarys')
    # plt.xlabel('Zeit [s]')
    # plt.ylabel('MEMS Verschiebung [nm]')
    # plt.legend()
    
    # plt.subplot(3,1,2)
    # plt.plot(P, M)
    # plt.plot(np.take(P, index), np.take(M, index), ls = '', marker = "o", label = 'segment boundarys')
    # plt.xlabel('Piezoposition [nm]')
    # plt.ylabel('MEMS Verschiebung [nm]')
    
    
    unload_Piezo, unload_MEMS = [],[]
    popt_log, pcov_log = [],[]
    S = []
    hmax = []
    for i in range(len(index_l)):
        up = P[index_h[i] : index_ul[i]]
        uM = M[index_h[i] : index_ul[i]]
        unload_Piezo.append(up[::-1])
        unload_MEMS.append(uM[::-1]) 
        hmax.append(np.max(uM))
        if i==len(index_l):
            par, cov = fitting(unload_Piezo[i] , unload_MEMS[i], [0.1, 0.95])
            #uncomment for power law fit
            #par, cov = fitting(unload_Piezo[i] , np.log(unload_MEMS[i]), [0.4, 0.95], (1.0, 1.0, 95), fit_func=func_exp)
        else:         
            par, cov = fitting(unload_Piezo[i] , unload_MEMS[i], [0.1, 0.95]) 
            #par, cov = fitting(unload_Piezo[i] , unload_MEMS[i], [0.3, 0.95], (1.0 ,1,0)) 
            #uncomment for power law fit
            #par, cov = fitting(unload_Piezo[i] , unload_MEMS[i], [0.4, 0.95], (1.0, 1.0 , 95), fit_func=func_exp)
        popt_log.append(par)
        pcov_log.append(cov)
        S.append(par[0])   
    plt.subplot(2,1,1)
    plt.plot(hmax, S, marker = '.', label='Messung'+str(j))
    plt.xlabel('maximale MEMS-Verschiebung [nm]')
    plt.ylabel('relative Steifigkeit')
    plt.legend()
    
    S_uload.append(S)
S_mean = []
S_std = []
S_uload=np.array(S_uload)

for i in range(len(S_uload[0])):
    S_mean.append(np.mean(S_uload[:,i]))
    S_std.append(np.std(S_uload[:,i]))

#plot stiffness versus indentation depth with errorbars 
plt.subplot(2,1,2)       
plt.errorbar(hmax, S_mean, yerr = S_std)
plt.xlabel('MEMS-Verschiebung [nm]')
plt.ylabel('Steifigkeit')
plt.title('S als Mittelwert mit 1 $\sigma$ Standarabweichung ')

    
    
