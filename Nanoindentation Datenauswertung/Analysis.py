# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 15:32:01 2022

@author: eim01
"""

from pre_processing_plain import *

fit_range = [0.6, 0.95]
#%%

start_index = 0 #Hertz fit range
end_index = 300
#PDMS Measurement

path = 'data/PDMS/EM-PDMS-50-10-50-0,7um'
Piezo, MEMS, time, Cap = imp_data(path)
Piezo_np_raw, MEMS_np_raw, time_np_raw, Cap_np_raw = data_conversion(Piezo, MEMS, time, Cap)

Piezo_np_raw = (Piezo_np_raw - Piezo_np_raw[0])*1000
Piezo_np, MEMS_np, time_np, Cap_np, poc_i = poc_detect(Piezo_np_raw, MEMS_np_raw, time_np_raw, Cap_np_raw)

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
Piezo, MEMS, time, Cap = imp_data(path)
Piezo_np_raw, MEMS_np_raw, time_np_raw, Cap_np_raw = data_conversion(Piezo, MEMS, time, Cap)

Piezo_np_raw = (Piezo_np_raw - Piezo_np_raw[0])*1000
Piezo_np, MEMS_np, time_np, Cap_np, poc_i = poc_detect(Piezo_np_raw, MEMS_np_raw, time_np_raw, Cap_np_raw)

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
plt.plot(np.take(Depth, index), np.take(Force, index), ls = '', marker = "o", label = 'segment boundarys')

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

path = ['data/PDMS/EM-PDMS-120-10-120-0,8um-1','data/PDMS/EM-PDMS-120-10-120-1,2um-1' ]

    
for i in path:
   Piezo, MEMS, time, Cap = imp_data(i) 
   Piezo_np, MEMS_np, time_np, Cap_np = data_conversion(Piezo, MEMS, time, Cap)
   Piezo_np = (Piezo_np - Piezo_np[0])*1000
   P, M, t, C, poc_i = poc_detect(Piezo_np, MEMS_np, time_np, Cap_np)
   
   Force = MEMS_np*K
   Depth = (Piezo_np-MEMS_np)
   Depth = Depth - Depth[poc_i]

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
#hier sind noch einige Fehler!!!
path = 'data/PDMS/EM-PDMS-120-10-120-1,2um-1'
Piezo, MEMS, time, Cap = imp_data(path)
Piezo_np, MEMS_np, time_np, Cap_np = data_conversion(Piezo, MEMS, time, Cap)
Piezo_np = (Piezo_np - Piezo_np[0])*1000

Force = MEMS_np*K
P, M, t, C, poc_i = poc_detect(Piezo_np, MEMS_np, time_np, Cap_np)
#poc_i = poc_i + find_nearest(Force[poc_i:np.argmax(Force)]) #point of contact moved to equilibrium force position F=0

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
popt_exp, pcov_exp = fitting(reversed_Depth, reversed_Force, fit_range, (0.1,1,0), fit_func=func_exp)
S = calc_stiff(popt_exp, reversed_Depth[-1])    #[nN/nm]
h_c = calc_hc(reversed_Depth[-1], reversed_Force[-1], S, eps=0.774) #[nm]
E_Op, E_reduced_Op = calc_emod(S, area_sphere(h_c))
idx_F0 = find_nearest(Force[poc_i:np.argmax(Force)])

E_hz_reduced, E_hz = calc_hertz(Force[poc_i:np.argmax(Force)], Depth[poc_i:np.argmax(Force)])



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
#plt.plot(Depth[poc_i:np.argmax(Force)], func_hertz(Depth[poc_i:np.argmax(Force)], E_hz_reduced), label='Hertz-Fit')
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

