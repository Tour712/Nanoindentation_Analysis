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

K = 3 # N/m
fit_range = [0.6, 0.95]
#%%

start_index = 0 #Hertz fit range
end_index = 300
#PDMS Measurement

path = 'data/PDMS/EM-PDMS-50-10-50-0,7um'
Piezo, MEMS, time = imp_data(path)
Piezo_np_raw, MEMS_np_raw, time_np_raw = data_conversion(Piezo, MEMS, time)

Piezo_np_raw = (Piezo_np_raw - Piezo_np_raw[0])*1000
Piezo_np, MEMS_np, time_np, poc_i = poc_detect(Piezo_np_raw, MEMS_np_raw, time_np_raw)

Depth_raw = (Piezo_np_raw-MEMS_np_raw)-Piezo_np_raw[poc_i]
Force_raw = MEMS_np_raw*K
Depth = (Piezo_np-MEMS_np)-Piezo_np[0]
Force = MEMS_np*K

#identify segment boundarys
#index_h = 948
#index_ul = 1441 
index_l, index_h,index_ul = data_splitting(Piezo_np, MEMS_np, time_np)   
unload_Depth = Depth[index_h : index_ul+1]
unload_Force = Force[index_h : index_ul+1] 

reversed_Depth = unload_Depth[::-1] #[nm]
reversed_Force = unload_Force[::-1] #[nN]

#calculations
popt_exp, pcov_exp = fitting(reversed_Depth, reversed_Force, fit_range, (0.1,1,0), fit_func=func_exp)
S = calc_stiff(popt_exp, reversed_Depth[-1])    #[nN/nm]
h_c = calc_hc(reversed_Depth[-1], reversed_Force[-1], S, eps=0.774) #[nm]
E_Op, E_reduced_Op = calc_emod(S, area_sphere(h_c))

E_hz_reduced, E_hz = calc_hertz(Force[start_index:end_index], Depth[start_index:end_index])

#####Plotting
plt.subplot(2,1,1)
plt.plot(Depth_raw, Force_raw, label='Messdaten')

plt.plot(reversed_Depth, func_exp(reversed_Depth, popt_exp[0], popt_exp[1], popt_exp[2]), label='power-law Fit')
plt.plot(Depth[start_index:end_index], func_hertz(Depth[start_index:end_index], E_hz_reduced), label='Hertz-Fit')
plt.xlabel('Tiefe [nm]')
plt.ylabel('Kraft [nN]')
plt.legend()

plt.subplot(2,1,2)
plt.plot(time_np_raw, Force_raw)
plt.xlabel('Zeit [s]')
plt.ylabel('Kraft [nN]')

#%%
start_index = 50 #Hertz fit range
end_index = 200
#PDMS Measurement

path = 'data/PDMS/EM-PDMS-50-10-50-0,5um-2'
Piezo, MEMS, time = imp_data(path)
Piezo_np_raw, MEMS_np_raw, time_np_raw = data_conversion(Piezo, MEMS, time)

Piezo_np_raw = (Piezo_np_raw - Piezo_np_raw[0])*1000
Piezo_np, MEMS_np, time_np, poc_i = poc_detect(Piezo_np_raw, MEMS_np_raw, time_np_raw)

Depth_raw = (Piezo_np_raw-MEMS_np_raw)-Piezo_np_raw[poc_i]
Force_raw = MEMS_np_raw*K
Depth = (Piezo_np-MEMS_np)-Piezo_np[0]
Force = MEMS_np*K

#identify segment boundarys
index_l, index_h,index_ul = data_splitting(Piezo_np, MEMS_np, time_np)  
index = np.array([index_l, index_h, index_ul]) 
unload_Depth = Depth[index_h : index_ul+1]
unload_Force = Force[index_h : index_ul+1] 

reversed_Depth = unload_Depth[::-1] #[nm]
reversed_Force = unload_Force[::-1] #[nN]

#calculations
popt_exp, pcov_exp = fitting(reversed_Depth, reversed_Force, fit_range, (0.1,1,0), fit_func=func_exp)
S = calc_stiff(popt_exp, reversed_Depth[-1])    #[nN/nm]
h_c = calc_hc(reversed_Depth[-1], reversed_Force[-1], S, eps=0.774) #[nm]
E_Op, E_reduced_Op = calc_emod(S, area_sphere(h_c))

E_hz_reduced, E_hz = calc_hertz(Force[start_index:end_index], Depth[start_index:end_index])

#####Plotting
plt.subplot(2,1,1)
plt.plot(Depth, Force, label='Messdaten')
#plt.plot(time_np ,MEMS_np, ls = '', marker = "+", markersize = 2)
plt.plot(np.take(Depth_raw, index), np.take(Force_raw, index), ls = '', marker = "o", label = 'segment boundarys')

plt.plot(reversed_Depth, func_exp(reversed_Depth, popt_exp[0], popt_exp[1], popt_exp[2]), label='power-law Fit')
plt.plot(Depth[start_index:end_index], func_hertz(Depth[start_index:end_index], E_hz_reduced), label='Hertz-Fit')
plt.xlabel('Tiefe [nm]')
plt.ylabel('Kraft [nN]')
plt.legend()

plt.subplot(2,1,2)
plt.plot(time_np_raw, Force_raw)
plt.xlabel('Zeit [s]')
plt.ylabel('Kraft [nN]')




#%%
#Plotting of multiple curves

path = []

for i in range(1,6):
    path.append('data/PDMS/EM-PDMS-50-10-50-0,5um'+'-'+str(i))
path.append('data/PDMS/EM-PDMS-50-10-50-0,3um')
#path.append('data/PDMS/EM-PDMS-50-10-50-0,8um')
path.append('data/PDMS/EM-PDMS-50-10-50-0,7um')
    
for i in path:
   Piezo, MEMS, time = imp_data(i) 
   Piezo_np_raw, MEMS_np_raw, time_np_raw = data_conversion(Piezo, MEMS, time)
   Piezo_np_raw = (Piezo_np_raw - Piezo_np_raw[0])*1000
   Piezo_np, MEMS_np, time_np, poc_i = poc_detect(Piezo_np_raw, MEMS_np_raw, time_np_raw)
   Depth_raw = (Piezo_np_raw-MEMS_np_raw)
   Force_raw = MEMS_np_raw*K
   plt.subplot(2,1,1)
   plt.plot(Depth_raw, Force_raw, label=i)
   plt.xlabel('Tiefe [nm]')
   plt.ylabel('Kraft [nN]')
   plt.legend()
   plt.subplot(2,1,2)
   plt.plot(time_np_raw, Force_raw, label=i)
   plt.xlabel('Zeit [s]')
   plt.ylabel('Kraft [nN]')
   plt.legend()


#%%

path = 'data/EM-SA-50-10-50-0,5um'
Piezo, MEMS, time = imp_data(path)
Piezo_np_raw, MEMS_np_raw, time_np_raw = data_conversion(Piezo, MEMS, time)
Piezo_np_raw = (Piezo_np_raw - Piezo_np_raw[0])*1000
Depth_raw = (Piezo_np_raw-MEMS_np_raw)-Piezo_np_raw[poc_i]
Force_raw = MEMS_np_raw*K
Piezo_np, MEMS_np, time_np, poc_i = poc_detect(Piezo_np_raw, MEMS_np_raw, time_np_raw)

Depth = (Piezo_np-MEMS_np)-Piezo_np_raw[poc_i]
Force = MEMS_np*K

plt.subplot(1,1,1)
plt.plot(Piezo_np_raw, MEMS_np_raw)
plt.xlabel('Piezopositzion [nm]')
plt.ylabel('MEMS Verschiebung [nm]')
plt.title('Messung auf Saphir-Probe')
#plt.legend()

# plt.subplot(2,1,2)
# plt.plot(Depth_raw, Force_raw)
# plt.xlabel('Depth [nm]')
# plt.ylabel('Kraft [nN]')
# plt.title('Messung auf Saphir-Probe')
# #plt.legend()