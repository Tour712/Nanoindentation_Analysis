# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 15:32:01 2022

@author: eim01
"""

from pre_processing_plain import *
# import numpy as np
# import csv
# import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit


def area_sphere(h_c, R=7500):
    A = np.pi*(2*R*h_c-h_c**2)
    return A

#h_c = np.linspace(0,100,100)
#plt.plot(h_c, area_sphere(7500, h_c))
#%%

start_index = 0
end_index = 300
#PDMS Measurement
K = 3 # N/m

path = 'data/EM-PDMS-50-10-50-0,7um'
Piezo, MEMS, time = imp_data(path)
Piezo_np, MEMS_np, time_np = data_conversion(Piezo, MEMS, time)

Piezo_np = (Piezo_np - Piezo_np[0])*1000
Piezo_np, MEMS_np, time_np, poc_i = poc_detect(Piezo_np, MEMS_np, time_np)

Depth = (Piezo_np-MEMS_np)-180
Force = MEMS_np*K

#identify segment boundarys
index_h = 918
index_ul = 1411    
unload_Depth = Depth[index_h : index_ul+1]
unload_Force = Force[index_h : index_ul+1] 

reversed_Depth = unload_Depth[::-1] #[nm]
reversed_Force = unload_Force[::-1] #[nN]

popt_exp, pcov_exp = fitting(reversed_Depth, reversed_Force, fit_range, (0.1,1,0), fit_func=func_exp)


S = calc_stiff(popt_exp, reversed_Depth[-1])    #[nN/nm]
print(S)
h_c = calc_hc(reversed_Depth[-1], reversed_Force[-1], S, eps=0.774) #[nm]
print(h_c)
E_Op, E_reduced_Op = calc_emod(S, area_sphere(h_c))


E_hertz = calc_hertz(Force[start_index:end_index], Depth[start_index:end_index])

plt.plot(Depth, Force)
plt.plot(reversed_Depth, func_exp(reversed_Depth, popt_exp[0], popt_exp[1], popt_exp[2]))
plt.plot(Depth[start_index:end_index], func_hertz(Depth[start_index:end_index], E_hertz))