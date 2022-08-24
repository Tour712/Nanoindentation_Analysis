# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 13:44:25 2022

@author: eim01
"""

from pre_processing_plain import *


path = 'data/MEMS Resolution Test/MEMS step test,0,5nm,10Hz-1'
Piezo, MEMS, time = imp_data(path)
Piezo_np_raw, MEMS_np_raw, time_np_raw = data_conversion(Piezo, MEMS, time)

Piezo_np_raw = (Piezo_np_raw - Piezo_np_raw[0])*1000


plt.plot(time_np_raw, Piezo_np_raw,ls = '-', markersize = 4, label='Piezo')
plt.plot(time_np_raw, MEMS_np_raw,ls = '-', markersize = 4, label='MEMS')
plt.xlabel('Zeit [s]')
plt.ylabel('Verschiebung von Piezo/MEMS [nm]')
plt.title('MEMS Stufentest')
plt.legend()

#%%
#evaluate drift data
path = 'data/Drift-2, pressure off, chuck not activated'
MEMS, time, T = imp_data(path)
MEMS_np_raw, time_np_raw, T_np_raw = data_conversion(MEMS, time, T)

#calculate 5min drift
m = 10  #drift Interval in min
drift = np.zeros(len(MEMS_np_raw))
for i in range(len(MEMS_np_raw)-m*2):
    drift[i] = (MEMS_np_raw[i]-MEMS_np_raw[i+m*2])/m
    
# plt.subplot(2,1,1)
# plt.plot(time_np_raw[15:], drift[15:])
# plt.xlabel('Zeit [min]')
# plt.ylabel('Driftrate [nm/min]')
# plt.title('Drift Test')
# plt.legend()

# plt.subplot(2,1,2)
# plt.plot(time_np_raw, MEMS_np_raw)
# plt.xlabel('Zeit [min]')
# plt.ylabel('Verschiebung von Piezo/MEMS [nm]')
# plt.title('Drift Test')
# plt.legend()

fig, ax1 = plt.subplots()
ax1.plot(time_np_raw, MEMS_np_raw)
ax1.set_xlabel('X-axis')
ax1.set_ylabel('Y1-axis')

ax2 = ax1.twinx()
ax2.plot(time_np_raw[0:1100], T_np_raw[0:1100], ls='-', markevery=100)
ax2.set_ylabel('Y2-axis')

plt.show()

#%%
#analysis of MEMS step response
path = ['data/MEMS Dynamik/1nm step','data/MEMS Dynamik/5nm step','data/MEMS Dynamik/10nm step','data/MEMS Dynamik/20nm step','data/MEMS Dynamik/50nm step', 'data/MEMS Dynamik/50nm step-2','data/MEMS Dynamik/100nm step','data/MEMS Dynamik/100nm step-2']
relax = []
for p in path:
    time, Piezo, MEMS = imp_data(p)
    time_np, Piezo_np, MEMS_np = data_conversion(time, Piezo, MEMS)
    relax.append(MEMS_np[5]-MEMS_np[750])
    
   # plt.plot(time_np[0:750], Piezo_np[0:750], label='Piezo')
    plt.plot(time_np[5:750], MEMS_np[5:750]-Piezo_np[5:750], label=p)
    plt.xlabel('Zeit [s]')
    plt.ylabel('MEMS Abweichung normalisiert [nm]')
    plt.yscale('linear')
    plt.title('Sprungantwort MEMS')
    plt.legend()


#plt.bar([1,5,10,20,50,50,100,100], relax)
#%%
#analysis of MEMS ramp response
path = 'data/MEMS Dynamik/1nm,100ms,0,4um'
time, Piezo, MEMS = imp_data(path)
time_np, Piezo_np, MEMS_np = data_conversion(time, Piezo, MEMS)
plt.plot(time_np, MEMS_np, label='MEMS')
plt.plot(time_np, Piezo_np, label='Piezo')
plt.xlabel('Zeit [s]')
plt.ylabel('MEMS/Piezo Verschiebung [nm]')
plt.yscale('linear')
plt.title('Rampe MEMS')
plt.legend()




