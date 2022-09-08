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