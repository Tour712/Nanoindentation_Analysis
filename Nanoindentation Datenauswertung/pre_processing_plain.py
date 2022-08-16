# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 16:03:06 2022

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

poc = 0     #point of contact
C_area = [24.5, 0]
e_tip = 0
nu_tip = 0

fit_range = [0.4, 0.95]  #


#############################################################################
#functions

def imp_data(path, column_names = ['Piezo_pos','Mems_displ','time']):
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
            #delete header lines           
    index = Piezo.index('eoh')
    del Piezo[0:index+1]   
    del MEMS[0:index+1]
    del time[0:index+1] 
    return Piezo, MEMS, time


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
       
    for index_l in range(index_start,len(Piezo_np)):
        piezo_del = Piezo_np[index_l]-Piezo_np[index_l-1]
        if (piezo_del) < 0.3:
            break

    for index_h in range(index_l,len(Piezo_np)):
        piezo_del = Piezo_np[index_h]-Piezo_np[index_h-1]
        if abs(piezo_del) > 0.3:
            break

    for index_ul in range(index_h,len(Piezo_np)):
        piezo_del = Piezo_np[index_ul]-Piezo_np[index_ul-1]
        if abs(piezo_del) < 0.3:
            break

    return index_l, index_h, index_ul
 

def area_func(h_c, C = [24.5, 0, 0]):
    A = C[0]*h_c**2 + C[1]*h_c**1 + C[2]*h_c**0.5
    return A        
   
#fitting functions
def func_lin(reversed_piezo, alpha, m ):
    return alpha*reversed_piezo+m

def func_exp(reversed_piezo, alpha, m, h_f):
    return alpha*(reversed_piezo - h_f)**m

def func_log(reversed_piezo, alpha, m, h_f ):
    return np.log(alpha)+m*np.log(reversed_piezo-h_f)
    
def func_hertz(Depth, E):
    return (4/3)*10**(-9)*E*(7500**0.5)*Depth**(3/2)

def fitting(reversed_piezo, reversed_MEMS, fit_range, *p0, fit_func=func_lin):
    start_index = int(len(reversed_piezo)*fit_range[0])
    end_index = int(len(reversed_piezo)*fit_range[1])
    return curve_fit(fit_func, reversed_piezo[start_index:end_index], reversed_MEMS[start_index:end_index], *p0, maxfev=10000)

def detect_poc(F, h, delta=3):
    for index, value in enumerate(F):
        if np.abs(value) > 0+delta:
            return index
        
def poc_detect(Piezo_np, MEMS_np, time_np):
    for index_poc, val in enumerate(MEMS_np):
        if (val) > 1.5:
            return Piezo_np[index_poc:], MEMS_np[index_poc:], time_np[index_poc:], index_poc
            break
        
def calc_stiff(curve_param, h_max):
    S = curve_param[1]*curve_param[0]*(h_max-curve_param[2])**(curve_param[1]-1)
    return S

def calc_hc(h_max, P_max, S, eps=0.762):
    h_s = eps*(P_max/S)
    h_c = h_max - h_s
    return h_c

def calc_emod(S, A, beta=1.05, nu_s=0.5, E_t=1140, nu_t=0.07):
    # S in [nN/nm], A in [nm^2]-> E in [nN/nm^2]
    # daher: E[nN/nm^2]*10^9= E[Pa]
    E_r = 10**9 *(S*np.sqrt(np.pi))/(2*beta*(np.sqrt(A)))
    E = (1- nu_s**2)/(1/E_r-(1-nu_t**2)/(E_t*10**9))
    return E, E_r

def calc_H(P_max, A):
    return (P_max*10**-3)/(A*10**-12)

#Hertz Analysis

def calc_hertz(P, h, nu_s=0.5, E_t=1140, nu_t=0.07):
    popt,pcov = curve_fit(func_hertz, h, P)
    E = (1- nu_s**2)/(1/popt-(1-nu_t**2)/(E_t*10**9))  
    return popt, E
