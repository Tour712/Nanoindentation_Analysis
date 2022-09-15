# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 10:24:23 2022

@author: eim01
"""

from functions import *
###########################################################################
#this script is for analysis of partial-unload measurement on Saphire

path = 'data/Saphir/PM-Sa-150ms, 3um, 20nms ,9cycles-3'
#path = 'data/Saphir/PM-Sa-150ms, 3um, 7nms ,3cycles-2'
#path = 'data/PDMS/12.09/PM-PDMS-150ms,3um,20nms,9cycles,1,5um offset'

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
        par, cov = fitting(unload_Piezo[i] , unload_MEMS[i], [0.5, 0.95])
        #uncomment for power law fit
        #par, cov = fitting(unload_Piezo[i] , np.log(unload_MEMS[i]), [0.4, 0.95], (1.0, 1.0, 95), fit_func=func_exp)
    else:         
        par, cov = fitting(unload_Piezo[i] , unload_MEMS[i], [0.5, 0.95]) 
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
    
    P, M, t, C = data_conversion(Piezo[j], MEMS[j], time[j], Cap[j])
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
 
    

#%%
###########################################################################
#this script is for analysis of partial-unload measurement on PDMS
path = 'data/PDMS/12.09/PM-PDMS-150ms,3um,20nms,9cycles,1,5um offset'

Piezo, MEMS, time, Cap = imp_data(path)
Piezo_np, MEMS_np, time_np, Cap_np = data_conversion(Piezo, MEMS, time, Cap)
Piezo_np = (Piezo_np - Piezo_np[0])*1000

Force = MEMS_np*K

P, M, t, C, poc_i = poc_detect(Piezo_np, MEMS_np, time_np, Cap_np)

Depth = Piezo_np - MEMS_np
Depth = Depth - Depth[poc_i]

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
    
plt.subplot(2,1,1)
plt.plot(time_np ,MEMS_np, ls = '', marker = "+", markersize = 2)
plt.plot(np.take(time_np, index), np.take(MEMS_np, index), ls = '', marker = "o", label = 'segment boundarys')
plt.xlabel('Zeit [s]')
plt.ylabel('MEMS Verschiebung [nm]')
plt.legend()

plt.subplot(2,1,2)
plt.plot(Piezo_np, MEMS_np)
plt.plot(np.take(Piezo_np, index), np.take(MEMS_np, index), ls = '', marker = "o", label = 'segment boundarys')
plt.xlabel('Piezoposition [nm]')
plt.ylabel('MEMS Verschiebung [nm]')



plt.figure()
plt.subplot(2,1,1)
plt.plot(Depth, Force)
plt.xlabel('Eindringtiefe')
plt.ylabel('Kraft')

plt.subplot(2,1,2)
#plt.plot(np.take(Depth, index), np.take(Force, index), ls = '', marker = "o", label = 'segment boundarys')
plt.xlabel('Eindringtiefe')
plt.ylabel('Kraft')

unload_Piezo, unload_MEMS = [],[]
popt_exp, pcov_exp = [],[]
popt_lin, pcov_lin = [],[]
popt_log, pcov_log = [],[]
S_log, S_exp = [], []
hmax = []
n=0
for i in range(len(index_l)-n):
    up = Depth[index_h[i+n] : index_ul[i+n]]
    uM = Force[index_h[i+n] : index_ul[i+n]]
    unload_Piezo.append(up[::-1])
    unload_MEMS.append(uM[::-1]) 
    hmax.append(np.max(up))
    plt.plot(up, uM)

    if i==(len(index_l)-1):
        par_l, cov_l = fitting(up[::-1] ,uM[::-1], [0.8, 0.95])
        par_log, cov_log = fitting(up[::-1] , np.log(uM[::-1]), [0.5, 1], fit_func=func_log)
        par_exp, cov_exp = fitting(up[::-1] ,uM[::-1] , [0.8, 0.95], fit_func=func_exp)
    else:         
        par_l, cov_l = fitting(up[::-1] ,uM[::-1], [0.1, 0.95])
        par_log, cov_log = fitting(up[::-1] , np.log(uM[::-1]), [0.1, 0.95], fit_func=func_log)
        par_exp, cov_exp = fitting(up[::-1] ,uM[::-1] , [0.1, 0.95], fit_func=func_exp)
        
        
    popt_exp.append(par_exp)
    popt_log.append(par_log)
    popt_lin.append(par_l[0])
    S_log.append(calc_stiff(par_log, up[0]))
    S_exp.append(calc_stiff(par_exp, up[0]))
    plt.plot(up, func_exp(up, par_log[0], par_log[1], par_log[2]),'r')
    

#%%
#This script is for the analysis of an array measurement on PDMS, with partial unload load function

path = 'data/PDMS/13.09/APM-PDMS-150ms,3um,20nms,9cycles,1,5um offset,5x1'
#path = 'data/PDMS/13.09/APM-PDMS-150ms,3um,20nms,9cycles,1,5um offset,5x1-3'
path = 'data/PDMS/13.09/APM-PDMS-150ms,3um,20nms,9cycles,1,5um offset,5x1-4'
path = 'data/PDMS/13.09/APM-PDMS-150ms,3um,20nms,9cycles,1,5um offset,9x1-5'
path = 'data/PDMS/14.09/APM-PDMS-150ms,3um,20nms,9cycles,1,5um offset,25%,9x1-6'
path = 'data/PDMS/14.09/APM-PDMS-150ms,4um,20nms,9cycles,1,5um offset,75%,5x1-7'
path = 'data/PDMS/14.09/APM-PDMS-150ms,4um,7nms,9cycles,1,5um offset,75%,5x1-8'
path = 'data/PDMS/14.09/APM-PDMS-150ms,4um,20nms,16cycles,1,5um offset,75%,16x1-9'

#path = 'data/PDMS/15.09/AR-PDMS-100ms,4um,20nms,1,5um offset,27x1-1'

Piezo_, MEMS_, time_, Cap_ = imp_data(path)
Piezo, MEMS, time, Cap, POC, X_val, Y_val = split_array(Piezo_, MEMS_, time_, Cap_)
S_load = []
S_uload = []
hmax = []
P_in = []
P_off = []
S_l, S_e, S_linear = [],[],[]
E_r_JKR, E_JKR = [],[]
f_err = 0
f_n = 0
s = 1 #skip first s measurement
for j in range(len(Piezo)-s):
    j=j+s
    
    P, M, t, C = data_conversion(Piezo[j], MEMS[j], time[j], Cap[j])
    Force = M*K
    P_, M_, t_, C_, poc = poc_detect(P, M, t, C)  
    P = (P-P[0] )*1000
    
    Depth = P-M
    Depth = Depth - Depth[poc]

    # store each measurement as np array in a list
    Piezo[j], MEMS[j], time[j] = P, M, t
       
    index, index_l, index_h, index_ul = [],[],[],[]
    
    index = data_splitting(P, M, t)
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
    
    # plt.figure()
    # plt.subplot(2,1,1)
    # plt.plot(t ,Force, ls = '', marker = "+", markersize = 2)
    # plt.plot(np.take(t, index), np.take(Force, index), ls = '', marker = "o", label = 'segment boundarys')
    # plt.xlabel('Zeit [s]')
    # plt.ylabel('Kraft [nN]')
    # plt.legend()
    
    # plt.subplot(2,1,2)
    # plt.plot(Depth, Force)
    # #plt.plot(np.take(Depth, index), np.take(Force, index), ls = '', marker = "o", label = 'segment boundarys')
    # plt.xlabel('Piezoposition [nm]')
    # plt.ylabel('Kraft [nN]')
    
    
    E_r_JKR.append(calc_JKRp(Depth, Force, R= 7500)[0])
    E_JKR.append(calc_JKRp(Depth, Force, R= 7500)[1])
    unload_Depth, unload_Force = [],[]
    hmax = []
    popt_exp, pcov_exp = [],[]
    popt_lin, pcov_lin = [],[]
    popt_log, pcov_log = [],[]
    S_log, S_exp = [], []
    hmax = []
    n=0
    for i in range(len(index_l)-n):
        up = Depth[index_h[i+n] : index_ul[i+n]]
        uM = Force[index_h[i+n] : index_ul[i+n]]
        unload_Depth.append(up[::-1])
        unload_Force.append(uM[::-1]) 
        hmax.append(np.max(up))
        #plt.plot(up, uM)
        f_n +=1
    
        if i==(len(index_l)-1):
            par_l, cov_l = fitting(up[::-1] ,uM[::-1], [0.85, 0.95])
            
            try:
                par_log, cov_log = fitting(up[::-1] , np.log(uM[::-1]), [0.85, 0.95], fit_func=func_log)
            except RuntimeError:
                print("Error - log_fit failed")
                continue
                
            try:
                par_exp, cov_exp = fitting(up[::-1] ,uM[::-1] , [0.85, 0.95], fit_func=func_exp)
            except RuntimeError:
                print("Error - power_fit failed")
                f_err +=1
                continue
        else:         
            par_l, cov_l = fitting(up[::-1] ,uM[::-1], [0.1, 0.95])
            try:
                par_log, cov_log = fitting(up[::-1] , np.log(uM[::-1]), [0.1, 0.95], fit_func=func_log)
            except RuntimeError:
                print("Error - log_fit failed")
                continue
            
            try:
                par_exp, cov_exp = fitting(up[::-1] ,uM[::-1] , [0.1, 0.95], fit_func=func_exp)
            except RuntimeError:
                print("Error - power_fit failed")
                f_err +=1
                continue
               
        
        #plt.plot(up[::-1], func_exp(up[::-1], par_exp[0], par_exp[1], par_exp[2]),'r')
        popt_exp.append(par_exp)
        popt_log.append(par_log)
        popt_lin.append(par_l[0])
        S_log.append(calc_stiff(par_log, up[0]))
        S_exp.append(calc_stiff(par_exp, up[0]))

    
    S_l.append(S_log)
    S_e.append(S_exp)
    S_linear.append(popt_lin)
print('number of failed fits:', f_err)
S_mean = [[],[],[]]
S_std = [[],[],[]]
cycles = np.arange(0,len(S_e[0]),1)
S_e = np.array(S_e)
S_l = np.array(S_l)
S_linear = np.array(S_linear)
for i in range(len(S_e[0])):
    S_mean[0].append(np.mean(S_e[:,i]))
    S_std[0].append(np.std(S_e[:,i]))
    
    S_mean[1].append(np.mean(S_l[:,i]))
    S_std[1].append(np.std(S_l[:,i]))
    
    S_mean[2].append(np.mean(S_linear[:,i]))
    S_std[2].append(np.std(S_linear[:,i]))

#plot stiffness versus indentation depth with errorbars 
plt.figure()       
plt.errorbar(cycles, S_mean[0], yerr = S_std[0], label='from power fit')
plt.errorbar(cycles, S_mean[1], yerr = S_std[1], label='from log fit')
plt.errorbar(cycles, S_mean[2], yerr = S_std[2], label='from linear fit')
plt.xlabel('MEMS-Verschiebung [nm]')
plt.ylabel('Steifigkeit')
plt.title('S als Mittelwert mit 1 $\sigma$ Standarabweichung ')
plt.legend()

