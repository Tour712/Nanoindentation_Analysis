# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 09:41:51 2022

@author: eim01
"""

from pre_processing_plain import *

def imp_data(path, column_names = ['Piezo_pos','Mems_displ','time', 'Cap']):
    '''
    read raw data into dictionary. Keywords are specified by the list 'column names'.
    Deletes Header lines
    
    return Parameters:
        Piezo, MEMS, time data as list
    '''
   #create empty list for storing measurement data
    Piezo = []
    MEMS = []
    time = []
    Cap = []
    with open(path) as f:
        raw_data = csv.DictReader(f, column_names, delimiter = '\t')    
        for line in raw_data:
            Piezo.append(line['Piezo_pos'])
            MEMS.append(line['Mems_displ'])
            time.append(line['time']) 
            Cap.append(line['Cap'])
            #delete header lines           
    index = Piezo.index('eoh')
    del Piezo[0:index+1]   
    del MEMS[0:index+1]
    del time[0:index+1] 
    del Cap[0:index+1]
    return Piezo, MEMS, time, Cap

def data_conversion(Piezo, MEMS, time, Cap):
    '''
    Parameters
    ----------
    Piezo : List of str
    MEMS : List of str
    time : List of str

    Returns
    -------
    Piezo: numpy array of float
    MEMS: numpy array of float
    time: numpy array of float

    '''
    for i in range(len(Piezo)):
        Piezo[i]= float(Piezo[i].replace(',','.'))   
        MEMS[i]= float(MEMS[i].replace(',','.')) 
        time[i]= float(time[i].replace(',','.'))
        Cap[i]= float(Cap[i].replace(',','.'))
    return np.array(Piezo), np.array(MEMS), np.array(time), np.array(Cap)

def calc_S(Cap_np_raw, Piezo_np_raw, index_l):
    S = (Cap_np_raw[index_l-1]-Cap_np_raw[0])/(Piezo_np_raw[index_l-1]-Piezo_np_raw[0])
    return S
    

path = 'data/S_calib/S_calib_1808-1'
Piezo, MEMS, time, Cap = imp_data(path)
Piezo_np_raw, MEMS_np_raw, time_np_raw, Cap_np_raw = data_conversion(Piezo, MEMS, time, Cap)

Piezo_np_raw = (Piezo_np_raw - Piezo_np_raw[0])*1000
Cap_np_raw = (Cap_np_raw - Cap_np_raw[0])*10**6

index_l, index_h,index_ul = data_splitting(Piezo_np_raw, MEMS_np_raw, time_np_raw) 
S = calc_S(Cap_np_raw, Piezo_np_raw, index_l)

#linear fit for load segment
popt_lin,cov = curve_fit(func_lin, Piezo_np_raw[0:index_l-1], Cap_np_raw[0:index_l-1], maxfev=10000)

#linear fit for unload segment
popt_lin2,cov2 = curve_fit(func_lin, Piezo_np_raw[index_h:][::-1], Cap_np_raw[index_h:][::-1], maxfev=10000)

# plt.plot(Piezo_np_raw, MEMS_np_raw)
# plt.plot([Piezo_np_raw[index_h]], [MEMS_np_raw[index_h]], marker = "o" )
plt.subplot(2,1,1)
plt.plot(Piezo_np_raw, Cap_np_raw, label ='Messdaten')
plt.plot(Piezo_np_raw, func_lin(Piezo_np_raw, popt_lin[0], popt_lin[1]), label ='Fit-Daten')
plt.xlabel('Piezoposition [nm]')
plt.ylabel('Kapazität [af]')
plt.grid()
plt.legend()

plt.subplot(2,1,2)
#plt.plot(Piezo_np_raw, Cap_np_raw, label ='Messdaten')
plt.plot(Piezo_np_raw[0:index_l-1], (Cap_np_raw[0:index_l-1] - func_lin(Piezo_np_raw[0:index_l-1], popt_lin[0], popt_lin[1]))/S)
plt.grid()
plt.xlabel('Piezoposition [nm]')
plt.ylabel('Residuen [nm]')
plt.legend()
 

 
