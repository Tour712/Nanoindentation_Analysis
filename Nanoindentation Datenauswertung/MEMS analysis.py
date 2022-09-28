# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 13:44:25 2022

@author: eim01
"""

from functions import *
import scipy.ndimage as ndi

plt.rcParams.update({'font.size': 14})
#Test step resolution
path = 'data/MEMS Resolution Test/MEMS step test,0,5nm,10Hz-1'
Piezo, MEMS, time, Cap = imp_data(path)
Piezo_np_raw, MEMS_np_raw, time_np_raw = data_conversion(Piezo, MEMS, time)

Piezo_np_raw = (Piezo_np_raw - Piezo_np_raw[0])*1000


plt.plot(time_np_raw, Piezo_np_raw, 'k', marker ='.', ls = '-', markersize = 4, label='Piezo')
plt.plot(time_np_raw, MEMS_np_raw, 'r', marker='x', ls = '-',alpha=0.8, markersize = 4, label='MEMS')
plt.xlabel('Zeit [s]')
plt.ylabel('Verschiebung von Piezo/MEMS [nm]')
plt.title('MEMS Stufentest')
plt.legend()

plt.figure()
plt.hist(MEMS_np_raw, 50)

#%%
#evaluate drift data
path = 'data/Drift-2, pressure off, chuck not activated'
MEMS, time, T, Cap = imp_data(path)
MEMS_np_raw, time_np_raw, T_np_raw = data_conversion(MEMS, time, T)
N = 100
T_rmean = ndi.uniform_filter1d(T_np_raw[0:1100], N, mode='constant', origin=-(N//2))[:-(N-1)]

#calculate 5min drift
m = 50  #drift Interval in min
drift_M = np.zeros(len(MEMS_np_raw))
for i in range(len(MEMS_np_raw)-m*2):
    drift_M[i] = (MEMS_np_raw[i]-MEMS_np_raw[i+m*2])/m
    
m = 10  #drift Interval in min
drift_T = np.zeros(len(T_np_raw))
for i in range(len(T_np_raw)-m*2):
    drift_T[i] = (T_np_raw[i]-T_np_raw[i+m*2])/m
    
plt.subplot(2,1,1)
plt.plot(time_np_raw[15:], drift_M[15:])
plt.xlabel('Zeit [min]')
plt.ylabel('Driftrate [nm/min]')
plt.title('Drift Test')
plt.legend()

plt.subplot(2,1,2)
plt.plot(time_np_raw, MEMS_np_raw)
plt.xlabel('Zeit [min]')
plt.ylabel('Verschiebung von Piezo/MEMS [nm]')
plt.title('Drift Test')
plt.legend()



fig, ax1 = plt.subplots()
ax1.plot(time_np_raw[0:len(T_rmean)], MEMS_np_raw[0:len(T_rmean)])
ax1.set_xlabel('Zeit [min]')
ax1.set_ylabel('MEMS-Verschiebung [nm]')

ax2 = ax1.twinx()
ax2.plot(time_np_raw[0:len(T_rmean)], T_rmean, 'r', ls='', marker='.')
ax2.set_ylabel('Temperatur [°C]')

plt.show()

#%%

#analysis of MEMS step response
path = ['data/MEMS Dynamik/10nm step']#,'data/MEMS Dynamik/1nm step','data/MEMS Dynamik/5nm step','data/MEMS Dynamik/20nm step','data/MEMS Dynamik/50nm step', 'data/MEMS Dynamik/50nm step-2','data/MEMS Dynamik/100nm step','data/MEMS Dynamik/100nm step-2']
relax = []
for p in path:
    time, Piezo, MEMS = imp_data(p)
    time_np, Piezo_np, MEMS_np = data_conversion(time, Piezo, MEMS)
    relax.append(MEMS_np[5]-MEMS_np[750])
    
    plt.plot(time_np[0:750], Piezo_np[0:750],linestyle='', marker='.', markersize=4, label='Piezo')
    plt.plot(time_np[0:750], MEMS_np[0:750],linestyle='', marker='.', markersize=4)#-Piezo_np[0:750], label=p)
    plt.xlabel('Zeit [s]')
    plt.ylabel('MEMS Abweichung normalisiert [nm]')
    plt.yscale('linear')
    plt.title('Sprungantwort MEMS')
    plt.legend()


#plt.bar([1,5,10,20,50,50,100,100], relax)
#%%
#analysis of MEMS resolution
path = 'data/MEMS Resolution Test/resolution in contact-6'
time, Piezo, MEMS = imp_data(path)
time_np, Piezo_np, MEMS_np = data_conversion(time, Piezo, MEMS)
MEMS_std = np.std(MEMS_np)
MEMS_mean = np.mean(MEMS_np)
Piezo_std = np.std(Piezo_np)
Piezo_mean = np.mean(Piezo_np)

plt.subplot(1,2,1)
plt.plot(time_np, MEMS_np,'r', marker='+',linestyle='-',linewidth =0.5, label='MEMS')
plt.plot(time_np, Piezo_np, 'k',marker='.',linestyle='-',linewidth =0.5, label='Piezo')
plt.ylim(-0.7,0.7)
#plt.text(0, 0.6, fontsize=14)
plt.xlabel('Zeit [s]')
plt.ylabel('MEMS/Piezo Verschiebung [nm]')
plt.yscale('linear')
#plt.title('file:'+path+'\n'+'MEMS std='+str(MEMS_std)+'    MEMS mean='+str(MEMS_mean)+'\n'+'Piezo std='+str(Piezo_std)+'   Piezo mean='+str(Piezo_mean))
plt.title('in Kontakt')
plt.legend()


path = 'data/MEMS Resolution Test/resolution out of contact-2'

time, Piezo, MEMS = imp_data(path)
time_np, Piezo_np, MEMS_np = data_conversion(time, Piezo, MEMS)
MEMS_std = np.std(MEMS_np)
MEMS_mean = np.mean(MEMS_np)
Piezo_std = np.std(Piezo_np)
Piezo_mean = np.mean(Piezo_np)

plt.subplot(1,2,2)
plt.plot(time_np, MEMS_np,'r', marker='+',linestyle='-',linewidth =0.5, label='MEMS')
plt.plot(time_np, Piezo_np, 'k',marker='.',linestyle='-',linewidth =0.5, label='Piezo')
plt.ylim(-0.7,0.7)
plt.xlabel('Zeit [s]')
#plt.ylabel('MEMS/Piezo Verschiebung [nm]')
plt.yscale('linear')
#plt.title('file:'+path+'\n'+'MEMS std='+str(MEMS_std)+'    MEMS mean='+str(MEMS_mean)+'\n'+'Piezo std='+str(Piezo_std)+'   Piezo mean='+str(Piezo_mean))
plt.title('in Luft')
#plt.legend()


#%%
#nur für Erzeugung der Abbildung
path = 'data/MEMS Resolution Test/resolution in contact-6'
time, Piezo, MEMS = imp_data(path)
time_np, Piezo_np, MEMS_np = data_conversion(time, Piezo, MEMS)
MEMS_std = np.std(MEMS_np)
MEMS_mean = np.mean(MEMS_np)
Piezo_std = np.std(Piezo_np)
Piezo_mean = np.mean(Piezo_np)

plt.subplot(2,2,1)
plt.plot(time_np, MEMS_np,'r', marker='x',linestyle='-',linewidth =0.5, label='MEMS')
plt.plot(time_np, Piezo_np, 'k',marker='.',linestyle='-',linewidth =0.5, label='Piezo')
plt.ylim(-0.7,0.7)
#plt.text(0, 0.6, fontsize=14)
plt.xlabel('Zeit [s]')
plt.ylabel('Verschiebung [nm]')
plt.yscale('linear')
#plt.title('file:'+path+'\n'+'MEMS std='+str(MEMS_std)+'    MEMS mean='+str(MEMS_mean)+'\n'+'Piezo std='+str(Piezo_std)+'   Piezo mean='+str(Piezo_mean))
plt.title('in Kontakt')
#plt.legend()


path = 'data/MEMS Resolution Test/resolution out of contact-2'

time, Piezo, MEMS = imp_data(path)
time_np, Piezo_np, MEMS_np = data_conversion(time, Piezo, MEMS)
MEMS_std = np.std(MEMS_np)
MEMS_mean = np.mean(MEMS_np)
Piezo_std = np.std(Piezo_np)
Piezo_mean = np.mean(Piezo_np)

plt.subplot(2,2,2)
plt.plot(time_np, MEMS_np,'r', marker='x',linestyle='-',linewidth =0.5, label='MEMS')
plt.plot(time_np, Piezo_np, 'k',marker='.',linestyle='-',linewidth =0.5, label='Piezo')
plt.ylim(-0.7,0.7)
plt.xlabel('Zeit [s]')
#plt.ylabel('MEMS/Piezo Verschiebung [nm]')
plt.yscale('linear')
#plt.title('file:'+path+'\n'+'MEMS std='+str(MEMS_std)+'    MEMS mean='+str(MEMS_mean)+'\n'+'Piezo std='+str(Piezo_std)+'   Piezo mean='+str(Piezo_mean))
plt.title('in Luft')
#plt.legend()



plt.rcParams.update({'font.size': 14})
#Test step resolution
path = 'data/MEMS Resolution Test/MEMS step test,0,5nm,10Hz-1'
Piezo, MEMS, time, Cap = imp_data(path)
Piezo_np_raw, MEMS_np_raw, time_np_raw = data_conversion(Piezo, MEMS, time)

Piezo_np_raw = (Piezo_np_raw - Piezo_np_raw[0])*1000

plt.subplot(2,1,2)
plt.plot(time_np_raw, Piezo_np_raw, 'k', marker ='.', ls = '-', markersize = 4, label='Piezo')
plt.plot(time_np_raw, MEMS_np_raw, 'r', marker='x', ls = '-',alpha=0.8, markersize = 4, label='MEMS')
plt.xlabel('Zeit [s]')
plt.ylabel('Verschiebung [nm]')
plt.title('MEMS Stufentest')
#plt.legend()

handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center')
