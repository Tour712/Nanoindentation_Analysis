# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 14:25:42 2022

@author: eim01
"""

from functions import *

def JKR_fit1(P, a_0, h_contact, P_adh):
    R = 7500
    return ((a_0**2)/R) * ((1+np.sqrt(1-P/P_adh)))**(4/3) -((2*a_0**2)/3*R) * ((1+np.sqrt(1-P/P_adh))/2)**(1/3) + h_contact
     
def DMT_fit(h, K, P_off):
    R = 7500
    return R**(0.5) *h**(3/2) *K -P_off
    
    
path = 'data/PDMS/EM-PDMS-120-10-120-1,2um-1'

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

#plot segment boundaries
index = [index_l, index_h, index_ul]

#popt_uload, pcov_uload = curve_fit(func_lin, reversed_Piezo, reversed_MEMS, (0.0,0.95))
#popt_uload, pcov_uload = curve_fit(JKR_fit1, reversed_Force, reversed_Depth, (100, 2200,-600), maxfev=20000, bounds = ([0,0, -700],[100,3000,-300]))
popt_DMT, pcov_DMT = curve_fit(DMT_fit, unload_Depth, unload_Force, (100,-600), maxfev=20000)

plt.subplot(2,1,1)
plt.plot(P, M, label = 'JKR')
plt.subplot(2,1,2)
plt.plot(Depth, Force, label = 'JKR')
plt.plot(reversed_Depth, reversed_Force)
#plt.plot(JKR_fit1(reversed_Force, popt_uload[0], popt_uload[1], popt_uload[2]), reversed_Force,  label = 'JKR')
plt.plot(reversed_Depth, DMT_fit(reversed_Depth, popt_DMT[0], popt_DMT[1]),  label = 'DMT')
plt.legend()

#%%


a=(np.linspace(0,100,1001))*10**-1 #
E_star = 3*10**6  #N/m^2
R = 10**(-5)   # m
P_off = - 1000 * 10**-9  # N
w_adh =  (-2/3) * (P_off/(np.pi*R)) 

h_jkr = ((a**2)/R - np.sqrt((2* np.pi *a * w_adh)/E_star))#*10**-6
F_jkr = ((4* E_star * a**3)/3*R - 2*np.sqrt(2* np.pi *E_star * w_adh*a**3)) #*10**-11

#a=(np.linspace(0,100,100))*10**-7  #
h_hertz = (a**2)/R
F_hertz =  (4/3)*E_star*(R**0.5)*h_hertz**(3/2)

plt.plot(h_jkr, F_jkr, label = 'JKR')
#plt.plot(h_hertz, F_hertz, label = 'Hertz; F in [N] und h in[m]')
plt.legend()


#%%
h = [-543.788, -534.501, -522.833, -505.644, -487.923, -474.96, -460.52, -442.36, -424.588, -408.419, -394.821, -377.307, -360.71, -346.718, -325.011, -300.86, -282.9, -262.037, -246.88, -224.863, -211.286, -193.259, -175.539, -154.332, -132.628, -113.743, -93.71, -72.768, -46.825, -21.77, 7.576, 30.174, 56.659, 84.041, 106.232, 128.223, 153.729, 183.145, 217.012, 246.052, 270.311, 299.043, 326.504, 345.172, 358.411]
P = [-0.144, -0.243, -0.354, -0.464, -0.589, -0.705, -0.818, -0.936, -1.052, -1.178, -1.284, -1.394, -1.501, -1.608, -1.743, -1.872, -1.972, -2.051, -2.112, -2.203, -2.273, -2.328, -2.406, -2.5, -2.58, -2.651, -2.741, -2.835, -2.944, -3.003, -3.041, -3.077, -3.13, -3.156, -3.163, -3.151, -3.105, -3.066, -2.98, -2.879, -2.8, -2.711, -2.586, -2.462, -2.336]

popt_jkr, cov = curve_fit(JKR_fit1, P, h, (100, 2200,-600), maxfev=20000, )
plt.plot(h, P)
plt.figure()
plt.plot(JKR_fit1(P, popt_jkr[0], popt_jkr[1], popt_jkr[2]))





