# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 15:32:01 2022

@author: eim01
"""

import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

############################################################################
   
def data_splitting(Piezo_np,MEMS_np,time_np):
    '''
    

    Parameters
    ----------
    Piezo_np : TYPE
        DESCRIPTION.
    MEMS_np : TYPE
        DESCRIPTION.
    time_np : TYPE
        DESCRIPTION.

    Returns
    -------
    index_l : int
        index of end of load segment.
    index_h : int
        index of end of hold segment.
    index_ul : int
        index of end of unload segment.

    '''
    
    for index_l in range(5,len(Piezo_np)):
        piezo_del = Piezo_np[index_l]-Piezo_np[index_l-1]
        if (piezo_del) < 0.3:
            break

    for index_h in range(index_l,len(Piezo_np)):
        piezo_del = Piezo_np[index_h]-Piezo_np[index_h-1]
        if abs(piezo_del) > 0.3:
            break

    for index_ul in range(index_h,len(Piezo_np)):
        piezo_del = Piezo_np[index_ul]-Piezo_np[index_ul-1]
        if abs(piezo_del) < 0.3:
            break
    
    return index_l, index_h, index_ul
 
def log_conversion(reversed_piezo, reversed_MEMS):
    pass
   
#fitting functions
def func_exp(reversed_piezo, alpha, m, h_f):
    return alpha*(reversed_piezo - h_f)**m

def func_log(reversed_piezo, alpha, m, h_f):
    return np.log(alpha)+m*np.log(reversed_piezo-h_f)
    
def func_lin(reversed_piezo, alpha, m, h_f):
    return alpha+m*(reversed_piezo-h_f)
    
    
def fitting(reversed_piezo, reversed_MEMS, fit_range, fit_func=func_exp):
    start_index = int(len(reversed_piezo)*(1-fit_range))
    end_index = len(reversed_piezo)
    return curve_fit(fit_func, reversed_piezo[start_index:end_index], reversed_MEMS[start_index:end_index])

############################################################################

fit_range=1
alpha_t = 0.1
m_t = 1.3
h_f = 20
noise = np.random.normal(0,1,101)

h_t = np.linspace(0, 100, 101)
P_t = alpha_t*(h_t-h_f)**m_t+noise


P_t_noise = P_t + noise

popt_te, pcov_te = fitting(h_t[50:], np.log(P_t[50:]), fit_range, fit_func=func_log)
popt_tl2, pcov_tl2 = fitting(h_t[50:], P_t[50:], fit_range, fit_func=func_lin)

# popt_te, pcov_te = curve_fit(func_expon, h_t, P_t, (1, 1, 5))
# popt_tl, pcov_tl = fitting(np.log(h_t), np.log(P_t), fit_range, fit_func=func_log)

#popt_t_noise, pcov_t_noise = curve_fit(func_exp, h_t, P_t_noise, (1, 1))
#print('Parameterfitting mit idealen Datensatz:', popt_t, '\nParameterfitting mit verrauschten Datensatz:', popt_t_noise)


plt.plot(h_t, P_t, label='Daten')
plt.plot(h_t, func_lin(h_t, popt_tl2[0], popt_tl2[1], popt_tl2[2]), label='lin fitting')
plt.plot(h_t, func_exp(h_t, popt_te[0], popt_te[1], popt_te[2]), label='log fitting')
plt.legend()
#plt.loglog()


