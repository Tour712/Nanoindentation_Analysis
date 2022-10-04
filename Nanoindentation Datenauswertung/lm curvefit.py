# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 09:41:29 2022

@author: eim01
"""
from functions import *
import numpy as np
from lmfit import Parameters,minimize, fit_report

# def JKR_fitting_lmfit(params,P,h):
#     a_0 = params['a_0']
#     h_contact = params['h_contact']
#     P_adh = params['P_adh']
#     R = params['R']
#     h_fit = ((a_0**2)/R) * ((1+np.sqrt(1-P/P_adh)))**(4/3) -((2*a_0**2)/3*R) * ((1+np.sqrt(1-P/P_adh))/2)**(1/3) + h_contact
#     return h_fit-h

def JKR_fitting_lmfit(params,P,h):
    a_R = params['a_R']
    h_contact = params['h_contact']
    P_adh = params['P_adh']
    h_fit = (a_R) * ((1+np.sqrt(1-P/P_adh)))**(4/3) -((2/3)*a_R) * ((1+np.sqrt(1-P/P_adh))/2)**(1/3) + h_contact
    return h_fit-h

# def JKR_fit1(P, a_0, h_contact, P_adh, R=7500):
#     return ((a_0**2)/R) * ((1+np.sqrt(1-P/P_adh)))**(4/3) -((2*a_0**2)/3*R) * ((1+np.sqrt(1-P/P_adh))/2)**(1/3) + h_contact
 
def JKR_fit1(P, a_R, h_contact, P_adh):
    return (a_R) * ((1+np.sqrt(1-P/P_adh)))**(4/3) -((2/3)*a_R) * ((1+np.sqrt(1-P/P_adh))/2)**(1/3) + h_contact
  
    
def DMT_fit(h, K, P_off):
    R = 7500
    return R**(0.5) *h**(3/2) *K -P_off

def calc_a0(a_R, R=7500):
    return np.sqrt(a_R*R)
    
def calc_delta(P_adh, R=7500):
    return -(P_adh*2)/(3*np.pi*R)

def calc_Er(delta, a_0, R=7500):
    E_r = (9*np.pi*R**2*delta)/(2*a_0**3)
    return E_r*10**3
    
path = 'C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS_10-1/test-5'
path = 'C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS_10-1/30.09/test-9'

Piezo_, MEMS_, time_, Cap_ = imp_data(path)
P, M, t, C = data_conversion(Piezo_, MEMS_, time_, Cap_)

P_, M_, t_, C_, poc = poc_detect(P, M, t, C)  
P = (P -P[0])*1000
Force = M * 3
Depth = (P - M)
#Depth = Depth - Depth[np.argmin(Force[0:np.argmax(Force)])]

index_l, index_h, index_ul = data_splitting(P, M, t)

load_Depth = Depth[poc:index_l]
load_Force = Force[poc:index_l]
unload_Depth = Depth[index_h : index_ul+1]
unload_Force = Force[index_h : index_ul+1] 

reversed_Depth = unload_Depth[::-1] #[nm]
reversed_Force = unload_Force[::-1] #[nN]

reversed_Depth = reversed_Depth[0:]
reversed_Force = reversed_Force[0:]

unload_Piezo = P[index_h : index_ul+1]
unload_MEMS = M[index_h : index_ul+1] 

reversed_Piezo = unload_Piezo[::-1] #[nm]
reversed_MEMS = unload_MEMS[::-1] #[nN]


# Defining the various parameters
params = Parameters()

params.add('a_R',value=2000, min=0,vary=True)
params.add('h_contact', value=2000,  vary=True)
params.add('P_adh', value=np.min(Force), max=np.min(Force), vary=True)
#params.add('R', value=1.0, min=0, vary=False)

# Calling the minimize function. Args contains the x and y data.
#fitted_params = minimize(linear_fitting_lmfit, params, args=(x,y,), method='least_squares')
fitted_params = minimize(JKR_fitting_lmfit, params, args=(reversed_Force, reversed_Depth), method='leastsq')

a_R = fitted_params.params['a_R'].value
h_c = fitted_params.params['h_contact'].value 
P_adh = fitted_params.params['P_adh'].value
#R = fitted_params.params['R'].value 

a_0 = calc_a0(a_R)
print(a_0)
delta = calc_delta(P_adh)
print(delta)
E_r = calc_Er(delta, a_0)
print(E_r)

E_rp, E_p = calc_JKRp(Depth, Force)
# Pretty printing all the statistical data
print(fit_report(fitted_params))

plt.subplot(3,1,1)
plt.plot(Depth, Force, label = 'JKR')
plt.plot(JKR_fit1(reversed_Force,a_R,h_c, P_adh), reversed_Force,  label = 'JKR fit')



plt.subplot(3,1,2)
plt.plot(reversed_Force, reversed_Depth, label = 'h(F)')
plt.plot(reversed_Force, JKR_fit1(reversed_Force,a_R,h_c, P_adh), label = 'h(F)')
plt.legend()

# a=np.linspace(0.0001,1,10)
# for i in a:
#     plt.figure()
#     plt.plot(JKR_fit1(reversed_Force,i,h_c, P_adh), reversed_Force,  label = 'JKR fit')