# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 09:41:29 2022

@author: eim01
"""
from functions import *
import numpy as np
from lmfit import Parameters,minimize, fit_report


# Define the fitting function
def linear_fitting_lmfit(params,x,y):
    m = params['m']
    c = params['c']
    y_fit = m*x + c
    return y_fit-y

def JKR_fitting_lmfit(params,P,h):
    a_0 = params['a_0']
    h_contact = params['h_contact']
    P_adh = params['P_adh']
    R = params['R']
    h_fit = ((a_0**2)/R) * ((1+np.sqrt(1-P/P_adh)))**(4/3) -((2*a_0**2)/3*R) * ((1+np.sqrt(1-P/P_adh))/2)**(1/3) + h_contact
    return h_fit-h

     
def DMT_fit(h, K, P_off):
    R = 7500
    return R**(0.5) *h**(3/2) *K -P_off
    
    
path = 'C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS_10-1/test-5'

Piezo_, MEMS_, time_, Cap_ = imp_data(path)
P, M, t, C = data_conversion(Piezo_, MEMS_, time_, Cap_)

P_, M_, t_, C_, poc = poc_detect(P, M, t, C)  
P = (P )*1000
Force = M * 3
Depth = (P - M)
Depth = Depth# - Depth[np.argmin(Force[0:np.argmax(Force)])]

index_l, index_h, index_ul = data_splitting(P, M, t)

unload_Depth = Depth[index_h : index_ul+1]
unload_Force = Force[index_h : index_ul+1] 

reversed_Depth = unload_Depth[::-1] #[nm]
reversed_Force = unload_Force[::-1] #[nN]

unload_Piezo = P[index_h : index_ul+1]
unload_MEMS = M[index_h : index_ul+1] 

reversed_Piezo = unload_Piezo[::-1] #[nm]
reversed_MEMS = unload_MEMS[::-1] #[nN]


# Defining the various parameters
params = Parameters()
# Slope is bounded between min value of 1.0 and max value of 3.0
params.add('a_0', min=0.0)
params.add('h_contact', value=3200)
params.add('P_adh', value=-2715, max=-2712)
params.add('R', value=7500, vary=False)

# Calling the minimize function. Args contains the x and y data.
#fitted_params = minimize(linear_fitting_lmfit, params, args=(x,y,), method='least_squares')
fitted_params = minimize(JKR_fitting_lmfit, params, args=(reversed_Force, reversed_Depth), method='least_squares')

a_0 = fitted_params.params['a_0'].value
h_c = fitted_params.params['h_contact'].value 
a_0 = fitted_params.params['a_0'].value
h_c = fitted_params.params['h_contact'].value 

# Getting the fitted values
#m = fitted_params.params['m'].value
#c = fitted_params.params['c'].value    

# Printing the fitted values
# print('The slope (m) is ', m)
# print('The intercept (c) is ', c)

# Pretty printing all the statistical data
print(fit_report(fitted_params))

plt.subplot(3,1,1)
plt.plot(Depth, Force, label = 'JKR')
plt.plot(reversed_Depth, reversed_Force)
#plt.plot(JKR_fit1(reversed_Force, popt_uload[0], popt_uload[1], popt_uload[2]), reversed_Force,  label = 'JKR')
#plt.plot(reversed_Depth, DMT_fit(reversed_Depth, popt_DMT[0], popt_DMT[1]),  label = 'DMT')


plt.subplot(3,1,2)
plt.plot(reversed_Force, reversed_Depth, label = 'h(F)')
plt.legend()

plt.subplot(3,1,3)
plt.plot(reversed_Force, reversed_Depth, label = 'h(F)')
plt.legend()