# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 09:41:51 2022

@author: eim01
"""


#from pre_processing_mod_mapping import split_array
from functions import *


    
#%%
# einfache Messung
plt.rc('xtick', labelsize=14) 
plt.rc('ytick', labelsize=14) 

path = 'data/S_calib/S_calib_0809-1'
#path = 'data/S_calib/S_calib_0509-1'
#path = 'data/S_calib/S_calib_0509-2'
path = 'data/S_calib/S_calib_0509-3'

Piezo, MEMS, time, Cap = imp_data(path)
Piezo_np_raw, MEMS_np_raw, time_np_raw, Cap_np_raw = data_conversion(Piezo, MEMS, time, Cap)

Piezo_np_raw = (Piezo_np_raw - Piezo_np_raw[0])*1000
Cap_np_raw = (Cap_np_raw )*10**6

index_l, index_h,index_ul = data_splitting(Piezo_np_raw, MEMS_np_raw, time_np_raw) 
S = calc_S(Cap_np_raw, Piezo_np_raw, index_l)

#linear fit for load segment
popt_lin,cov = curve_fit(func_lin, Piezo_np_raw[0:index_l-1], Cap_np_raw[0:index_l-1], maxfev=10000)

#linear fit for unload segment
popt_lin2,cov2 = curve_fit(func_lin, Piezo_np_raw[index_h:][::-1], Cap_np_raw[index_h:][::-1], maxfev=10000)


fig,ax1 = plt.subplots()
ax3 = ax1.twinx()

ax1.plot(Piezo_np_raw, Cap_np_raw*10**-6,'k', label ='Messdaten')
#ax1.plot(Piezo_np_raw, func_lin(Piezo_np_raw, popt_lin[0], popt_lin[1])*10**-6, label ='Fit-Daten')
ax1.set_xlabel('Piezoposition [nm]',fontsize=14)
ax1.set_ylabel('Kapazität [pF]', fontsize=14)
ax1.grid()
#ax1.legend()

ax3.plot(Piezo_np_raw[0:index_l][::10], ((Cap_np_raw[0:index_l] - func_lin(Piezo_np_raw[0:index_l], popt_lin[0], popt_lin[1]))/S)[::10],'r', marker='x', label ='Residuen Ladesegment')
ax3.plot(Piezo_np_raw[index_h:index_ul][::10], ((Cap_np_raw[index_h:index_ul] - func_lin(Piezo_np_raw[index_h:index_ul], popt_lin[0], popt_lin[1]))/S)[::10],'g', marker='*', label ='Residuen Entladesegment')
#plot an empty array with same style as data, to show all labels in one legend
ax3.plot([],[],'k', label ='Messdaten')

ax3.set_ylabel('Residuen [nm]',fontsize=14)
ax3.legend()


#plt.plot(Piezo_np_raw, Cap_np_raw, label ='Messdaten')
# ax2.plot(Piezo_np_raw[0:index_l-1], (Cap_np_raw[0:index_l-1] - func_lin(Piezo_np_raw[0:index_l-1], popt_lin[0], popt_lin[1]))/S)
# ax2.grid()
# ax2.set_xlabel('Piezoposition [nm]',fontsize=14)
# ax2.set_ylabel('Residuen [nm]',fontsize=14)
# ax2.legend()


#%%
# Arraymessung

path = 'data/S_calib/S_calib_2308-3x3'
#path = 'data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-3'
Piezo_, MEMS_, time_, Cap_ = imp_data(path)
Piezo, MEMS, time, Cap, POC, X_val, Y_val = split_array(Piezo_, MEMS_, time_, Cap_)

Results = []
indent_depth = []
fit_param = []

S = []
S_fit_load = []
S_fit_uload = []
for i, e in enumerate(Piezo):
    
    P, M, t, C = data_conversion(Piezo[i], MEMS[i], time[i], Cap[i])
    # Piezo offset position and conversion to [nm]
    P = (P - P[0])*1000
    C = (C - C[0])*10**6
    
    # detect POC
    P_r, M_r, t_r, C_r,p = poc_detect(P, M, t, C)  
    # store each measurement as np array in a list
    Piezo[i], MEMS[i], time[i], Cap[i] = P_r, M_r, t_r, C_r
    
    index_l, index_h,index_ul = data_splitting(P_r, M_r, t_r) 
    S.append(calc_S(C_r, P_r, index_l)) 
    plt.plot(P_r, C_r)
    
    #linear fit for load segment
    popt_lin,cov = curve_fit(func_lin, P_r[0:index_l-1], C_r[0:index_l-1], maxfev=10000)
    #linear fit for unload segment
    popt_lin2,cov2 = curve_fit(func_lin, P_r[index_h:][::-1], C_r[index_h:][::-1], maxfev=10000)
    S_fit_load.append(popt_lin[0])
    S_fit_uload.append(popt_lin2[0])
    plt.plot()
