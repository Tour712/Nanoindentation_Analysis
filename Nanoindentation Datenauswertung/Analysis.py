# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 15:32:01 2022

@author: eim01
"""

from functions import *

fit_range = [0.6, 0.95]

#%%
#Plotting of multiple curves

path = ['data/PDMS/EM-PDMS-120-10-120-0,8um-1','data/PDMS/EM-PDMS-120-10-120-1,2um-1' ]

path = ['data/PDMS/12.09/EM-PDMS-110-10-110, 20nms,2,2um-2', 'data/PDMS/12.09/EM-PDMS-125-10-125, 20nms,2,5um-3', 'data/PDMS/12.09/EM-PDMS-150-10-150, 20nms,3um-4']#,'data/PDMS/13.09/EM-PDMS-150-10-150-20nms-3um-5']

    
for i in path:
   Piezo, MEMS, time, Cap = imp_data(i) 
   Piezo_np, MEMS_np, time_np, Cap_np = data_conversion(Piezo, MEMS, time, Cap)
   Piezo_np = (Piezo_np - Piezo_np[0])*1000
   P, M, t, C, poc_i = poc_detect(Piezo_np, MEMS_np, time_np, Cap_np)
   Force = MEMS_np*K
   Depth = (Piezo_np-MEMS_np)
   Depth = Depth - Depth[poc_i]
   print(calc_JKRp(Depth, Force, R= 7500))
   plt.subplot(2,1,1)
   plt.plot(Piezo_np, Force, label=i)
   plt.xlabel('Piezo [nm]')
   plt.ylabel('Kraft [nN]')
   plt.grid(b=True)
   plt.legend()
   plt.subplot(2,1,2)
   plt.plot(Depth, Force, label=i)
   plt.xlabel('Tiefe [nm]')
   plt.ylabel('Kraft [nN]')
   plt.grid(b=True)
   plt.legend()


#%%
#Analysis of PDMS measurement (einfache Messung)
#hier sind noch einige Fehler!!!
#path = 'data/PDMS/EM-PDMS-120-10-120-1,2um-1'
path = 'data/PDMS/09.09/EM-PDMS-110-10-110, 20nms,2,2um'
path = 'data/PDMS/12.09/EM-PDMS-110-10-110, 20nms,2,2um-2'
path = 'data/PDMS/12.09/EM-PDMS-125-10-125, 20nms,2,5um-3'
path = 'data/PDMS/13.09/EM-PDMS-150-10-150-20nms-3um-5'
#path = 'data/PDMS/12.09/EM-PDMS-150-10-150, 20nms,3um-4'
Piezo, MEMS, time, Cap = imp_data(path)
Piezo_np, MEMS_np, time_np, Cap_np = data_conversion(Piezo, MEMS, time, Cap)
Piezo_np = (Piezo_np - Piezo_np[0])*1000

Force = MEMS_np*K
P, M, t, C, poc_i = poc_detect(Piezo_np, MEMS_np, time_np, Cap_np)
#poc_i = poc_i + find_nearest(Force[poc_i:np.argmax(Force)]) #point of contact moved to equilibrium force position F=0, necessary for Hertz analysis

Depth = Piezo_np - MEMS_np
Depth = Depth - Depth[poc_i]

E_r_jkr, E_jkr = calc_JKRp(Depth, Force, R= 7500)
print(E_r_jkr)

#identify segment boundarys
index_l, index_h,index_ul = data_splitting(Piezo_np, MEMS_np, time_np)  
index = np.array([poc_i, index_l, index_h, index_ul]) 
unload_Depth = Depth[index_h : index_ul+1]
unload_Force = Force[index_h : index_ul+1] 

reversed_Depth = unload_Depth[::-1] #[nm]
reversed_Force = unload_Force[::-1] #[nN]

# #calculations
popt_exp, pcov_exp = fitting(reversed_Depth, reversed_Force, fit_range, fit_func=func_exp)
S = calc_stiff(popt_exp, reversed_Depth[-1])    #[nN/nm]
h_c = calc_hc(reversed_Depth[-1], reversed_Force[-1], S, eps=0.774) #[nm]
E_Op, E_reduced_Op = calc_emod(S, area_sphere(h_c))
idx_F0 = find_nearest(Force[poc_i:np.argmax(Force)])
#E_hz_reduced, E_hz  = calc_hertz(Force[poc_i:index_l], Depth[poc_i:index_l])



plt.subplot(2,1,1)
plt.plot(Piezo_np, MEMS_np)
plt.xlabel('Piezopositzion [nm]')
plt.ylabel('MEMS Verschiebung [nm]')
plt.grid()
plt.title(path)


plt.subplot(2,1,2)
plt.plot(Depth, Force)
plt.plot(np.take(Depth, index), np.take(Force, index), ls = '', marker = "o", label = 'segment boundarys')
plt.plot(reversed_Depth, func_exp(reversed_Depth, popt_exp[0], popt_exp[1], popt_exp[2]), label='power-law Fit')
#plt.plot(Depth[poc_i:index_l], func_hertz(Depth[poc_i:index_l], E_hz_reduced), label='Hertz-Fit')
plt.xlabel('Depth [nm]')
plt.ylabel('Kraft [nN]')
plt.grid()
plt.legend()



#%%
# tip-approach graph
plt.figure()
plt.plot(Piezo_np[0:index_l], MEMS_np[0:index_l]*3,'r', marker='.',markersize=2, linestyle='')
plt.xlabel('Piezopositzion [nm]', size=18)
plt.ylabel('Kraft [nN]', size=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid()

