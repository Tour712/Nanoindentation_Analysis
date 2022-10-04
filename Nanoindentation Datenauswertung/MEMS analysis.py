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
from matplotlib.ticker import FormatStrFormatter
path = 'C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/Drift/Drift-1, pressure on, chuck not activated'
path = 'data/Drift-2, pressure off, chuck not activated'
#path = 'data/Drift/Drift-3, vacuum chuck activated, on 3inchwafer'
MEMS, time, T, Cap = imp_data(path)
MEMS_np_raw, time_np_raw, T_np_raw = data_conversion(MEMS, time, T)
s = 5 #start at s min 
e = 571 #end at e min
MEMS_np_raw, time_np_raw, T_np_raw = MEMS_np_raw[s*2:e*2], time_np_raw[s*2:e*2], T_np_raw[s*2:e*2]
time_np_raw =time_np_raw-time_np_raw[0]
N = 30
T_rmean = ndi.uniform_filter1d(T_np_raw, N, mode='constant', origin=-(N//2))[:-(N-1)]

#calculate 5min drift
m = 10  #drift Interval in min
drift_M = np.zeros(len(MEMS_np_raw))
for i in range(len(MEMS_np_raw)-m*2):
    drift_M[i] = (MEMS_np_raw[i]-MEMS_np_raw[i+m*2])/m
    
m = 10  #drift Interval in min
drift_T = np.zeros(len(T_np_raw))
for i in range(len(T_rmean)-m*2):
    drift_T[i] = (T_rmean[i]-T_rmean[i+m*2])/m
    
# plt.subplot(2,1,1)
# plt.plot(time_np_raw, drift_T)
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



fig, (ax1,ax3) = plt.subplots(2,1)
ax1.plot(time_np_raw[0:len(T_rmean)], MEMS_np_raw[0:len(T_rmean)],'k', linestyle='', marker='.')
ax1.set_xlabel('Zeit [min]')
ax1.set_ylim(124,146)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
ax1.set_ylabel('Verschiebung [nm]')
ax1.set_xlim(-5,500)
ax1.grid(visible=True)
ax2 = ax1.twinx()
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
ax2.plot(time_np_raw[0:len(T_rmean)], T_rmean, 'r', ls='', marker='1',label='Temperatur')
ax2.set_ylabel('Temperatur [°C]')
ax2.plot([],'k',linestyle='', marker='.', label='Verschiebung')
ax2.legend(loc='lower right',prop={'size': 10})

ax3.plot(time_np_raw, drift_M, 'k', label ='Positionsdrift')
ax3.grid(visible=True)
ax3.set_xlabel('Zeit [min]')
ax3.set_ylabel('Positionsdrift [nm/min]')
ax3.set_ylim(-0.15,0.15)
ax3.set_xlim(-5,500)
ax4 = ax3.twinx()
ax4.plot(time_np_raw, drift_T*1000, 'r', label='Temperaturdrift')
ax4.set_ylabel('Temperaturdrift [mK/min]')
ax4.plot([],'k', label='Positionsdrift')
ax4.set_ylim(-0.17,0.17)
ax4.legend(prop={'size': 10})
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


#%%
#Piezotisch Dynamik
path = 'data/Piezotisch Dynamik/step response 10nm.tmp.txt'
path = ['data/Piezotisch Dynamik/step response 1nm.tmp.txt','data/Piezotisch Dynamik/step response 10nm.tmp.txt', 
        'data/Piezotisch Dynamik/step response 30nm.tmp.txt','data/Piezotisch Dynamik/step response 50nm.tmp.txt',
        'data/Piezotisch Dynamik/step response 70nm.tmp.txt','data/Piezotisch Dynamik/step response 100nm.tmp.txt']
plt.rcParams.update({'font.size': 14})
def func_pt1(t, A, tau,y_0):
    return A*(1-np.exp(-t/tau))+y_0
popt,cov = [], []
colors=['g','b','r','k','c','y','w']
for a,j in enumerate(path):
    P = []
    t = []
    
    with open(j) as f:
        raw_data = csv.DictReader(f, ['time', 'Piezo_pos'], delimiter = '\t')    
        for line in raw_data:
            P.append(line['Piezo_pos'])
            t.append(line['time'])          
    P, t = P[1:], t[1:]
    for i in range(len(P)):
        P[i]= float(P[i].replace(',','.'))   
        t[i]= float(t[i].replace(',','.'))
    
    
    P, t = np.array(P), np.array(t)
    P = (P - P[0])*1000
    popt.append(curve_fit(func_pt1, t, P)[0])
    plt.subplot(1,2,1)
    plt.plot(t, P, c=colors[a])  
    plt.plot(t, func_pt1(t, *popt[a]),linestyle='dashed',c='gray')
    plt.grid(visible=True)
    #plt.yscale('log')
    plt.xlim([0,20])
    plt.xlabel('Zeit [ms]')
    plt.ylabel('Piezo [nm]')
    
    plt.subplot(1,2,2)
    plt.plot(t, P, c=colors[a])   
    plt.grid(visible=True)
    plt.yscale('log')
    plt.xlim([0,20])
    plt.xlabel('Zeit [ms]')
    plt.ylabel('Piezo [nm]')
    
popt = np.array(popt)
print(popt[:,1])
plt.subplot(1,2,1)
plt.plot([],[], linestyle='dashed',c='gray', label='Fit')
# plt.plot(t, func_pt1(t, *popt[2]),linestyle='dashed')
plt.legend()
