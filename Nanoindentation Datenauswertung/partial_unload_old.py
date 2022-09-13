# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 15:15:32 2022

@author: eim01
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 16:03:06 2022

@author: eim01
"""
from functions import *

#create empty list for storing measurement data
Piezo = []
MEMS = []
time = []

poc = 0     #point of contact
area_func = []
e_tip = 0
nu_tip = 0

fit_range = [0, 1]  #
#############################################################################
#partial unload
#data import and conversion
path = 'data/Test Measurement/partial unload on sapphire-4'
Piezo, MEMS, time, Cap = imp_data(path)
Piezo_np, MEMS_np, time_np, Cap_np = data_conversion(Piezo, MEMS, time, Cap)

#Piezo offset position and conversion to [nm]
Piezo_np = (Piezo_np - Piezo_np[0])*1000
Piezo_np, MEMS_np, time_np, poc_i = poc_detect(Piezo_np, MEMS_np, time_np)
Data = np.array([Piezo_np, MEMS_np, time_np])   #3xn array containing all the data

#############################################################################
#split loading curve in load-, hold- and unload segment

#identify segment boundarys
index, index_l, index_h, index_ul = [],[],[],[]

index = data_splitting(Piezo_np, MEMS_np, time_np)
index_l.append(data_splitting(Piezo_np, MEMS_np, time_np)[0])
index_h.append(data_splitting(Piezo_np, MEMS_np, time_np)[1])
index_ul.append(data_splitting(Piezo_np, MEMS_np, time_np)[2]) 


while index[-1] < np.argmin(MEMS_np):
    ind = []
    ind = data_splitting(Piezo_np, MEMS_np, time_np, index_start= index[-1]+1)
    index = index + ind
    index_l.append(ind[0])
    index_h.append(ind[1])
    index_ul.append(ind[2])
# index = [x - 1 for x in index]
# index_l = [x - 1 for x in index_l]
# index_h = [x - 1 for x in index_h]
# index_ul = [x - 1 for x in index_ul]

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

fig, ax = plt.subplots()
ax.plot(Piezo_np, MEMS_np)
ax.scatter(np.take(Piezo_np, index), np.take(MEMS_np, index))

for i, txt in enumerate(index):
    ax.annotate(i, (np.take(Piezo_np, index)[i], np.take(MEMS_np, index)[i]))


#lists of numpy arrays, containing the unloading segments, reverse order of data points
#fit curve to every segment and store fitting results in list
unload_Piezo, unload_MEMS = [],[]
popt_log, pcov_log = [],[]
S = []
for i in range(len(index_l)):
    up = Piezo_np[index_h[i] : index_ul[i]]
    #up =  up[::-1]
    uM = MEMS_np[index_h[i] : index_ul[i]]
    #uM =  uM[::-1]
    unload_Piezo.append(up[::-1])
    unload_MEMS.append(uM[::-1]) 
    if i==len(index_l):
        par, cov = fitting(unload_Piezo[i] , np.log(unload_MEMS[i]), [0.3, 0.95], (1.0,1,0), fit_func=func_log)
        #uncomment for power law fit
        #par, cov = fitting(unload_Piezo[i] , np.log(unload_MEMS[i]), [0.4, 0.95], (1.0, 1.0, 95), fit_func=func_exp)
    else:         
        par, cov = fitting(unload_Piezo[i] , np.log(unload_MEMS[i]), [0.3, 0.95], (1.0 ,1,0), fit_func=func_log) 
        #uncomment for power law fit
        #par, cov = fitting(unload_Piezo[i] , unload_MEMS[i], [0.4, 0.95], (1.0, 1.0 , 95), fit_func=func_exp)
    popt_log.append(par)
    pcov_log.append(cov)
    ax.plot(unload_Piezo[i], func_exp(unload_Piezo[i], par[0], par[1], par[2]), label = 'log fit' + str(i))
    S.append(calc_stiff(par, up[0]))


        
# #############################################################################
# #curve fitting

# popt_exp, pcov_exp = fitting(reversed_piezo, reversed_MEMS, fit_range, (0.1,1,0), fit_func=func_exp)
# popt_lin, pcov_lin = fitting(reversed_piezo, reversed_MEMS, fit_range, fit_func=func_lin)


# #plot fitting result and raw data
# plt.plot(reversed_piezo, reversed_MEMS, label='data')
# plt.plot(reversed_piezo, func_exp(reversed_piezo, popt_log[0], popt_log[1], popt_log[2]), label = 'log fit')
# plt.plot(reversed_piezo, func_exp(reversed_piezo, popt_exp[0], popt_exp[1], popt_exp[2]), label = 'power-law fit')
# plt.legend()

# #############################################################################




