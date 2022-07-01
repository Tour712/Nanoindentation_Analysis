# -*- coding: utf-8 -*-
"""
Created on Thu Jun  29 13:00:06 2022

@author: eim01
"""


import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#create empty list for storing measurement data
Piezo = []
MEMS = []
time = []
Cap = []

poc = 0     #point of contact
area_func = []
e_tip = 0
nu_tip = 0

fit_range = [0, 1]  #


#%%
#############################################################################
#functions

def imp_data(path, column_names = ['Piezo_pos','Mems_displ','time','Cap']):
    '''
    read raw data into dictionary. Keywords are specified by the list 'column names'.
    Deletes Header lines
    
    return Parameters:
        Piezo, MEMS, time data as list
    '''
   
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


def data_conversion(Piezo, MEMS, time):
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
    return np.array(Piezo), np.array(MEMS), np.array(time)

def poc_detect(Piezo_np, MEMS_np, time_np):
    for index_poc, val in enumerate(MEMS_np):
        if (val) > 1.5:
            return Piezo_np[index_poc:], MEMS_np[index_poc:], time_np[index_poc:], index_poc
            break
    
    
def data_splitting(Piezo_np,MEMS_np,time_np, index_start = 5):
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

    for index_l in range(index_start, len(Piezo_np)):
        piezo_del = Piezo_np[index_l]-Piezo_np[index_l-1]
        if (piezo_del) < 0.3:
            #index_load.append(index_l)
            break

    for index_h in range(index_l, len(Piezo_np)):
        piezo_del = Piezo_np[index_h]-Piezo_np[index_h-1]
        if abs(piezo_del) > 0.3:
            #index_hold.append(index_h)
            break

    for index_ul in range(index_h, len(Piezo_np)):
        piezo_del = Piezo_np[index_ul]-Piezo_np[index_ul-1]
        if abs(piezo_del) < 0.3:
            break
    return [index_l, index_h, index_ul]           
    
   
#fitting functions
def func_lin(reversed_piezo, alpha, m ):
    return alpha*reversed_piezo+m

def func_exp(reversed_piezo, alpha, m, h_f):
    return alpha*(reversed_piezo - h_f)**m

def func_log(reversed_piezo, alpha, m, h_f ):
    return np.log(alpha)+m*np.log(reversed_piezo-h_f)
    
def fitting(reversed_piezo, reversed_MEMS, fit_range, *p0, fit_func=func_lin):
    start_index = int(len(reversed_piezo)*fit_range[0])
    end_index = int(len(reversed_piezo)*fit_range[1])
    return curve_fit(fit_func, reversed_piezo[start_index:end_index], reversed_MEMS[start_index:end_index], *p0, maxfev=5000 )
    
def calc_stiff(curve_param, h_max):
    S = curve_param[1]*curve_param[0]*(h_max-curve_param[2])**(curve_param[1]-1)
    return S

def calc_emod(S, A, beta):
    E = S*np.sqrt(np.pi)/2*beta*(np.sqrt(A))
    return E

def calc_hf(reversed_piezo, reversed_MEMS,fit_range_hf):
    pass
    

#############################################################################
#partial unload
#data import and conversion
path = 'data/modulus mapping on saphir, 2x2'
Piezo, MEMS, time, Cap = imp_data(path)
print(MEMS[0:15])
#Piezo_np, MEMS_np, time_np = data_conversion(Piezo, MEMS, time)

#Piezo offset position and conversion to [nm]
# Piezo_np = (Piezo_np - Piezo_np[0])*1000
# Piezo_np, MEMS_np, time_np, poc_i = poc_detect(Piezo_np, MEMS_np, time_np)
# Data = np.array([Piezo_np, MEMS_np, time_np])   #3xn array containing all the data

# #############################################################################
# #split loading curve in load-, hold- and unload segment

# #identify segment boundarys
# index, index_l, index_h, index_ul = [],[],[],[]



