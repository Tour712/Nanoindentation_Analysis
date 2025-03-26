# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 16:03:06 2022

@author: eim01
"""

import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lmfit import Parameters,minimize, fit_report

plt.rcParams.update({'font.size': 14})
poc = 0     #point of contact
C_area = [24.5, 0]
e_tip = 0
nu_tip = 0
K = 3.01 # N/m
fit_range = [0.4, 0.95]  #


#############################################################################
#functions

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
    
    if Cap[10] is None:
        return Piezo, MEMS, time
    else:
        return Piezo, MEMS, time, Cap


def data_conversion(Piezo, MEMS, time, Cap=0):
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
    if Cap == 0:
        for i in range(len(Piezo)):
            Piezo[i]= float(Piezo[i].replace(',','.'))   
            MEMS[i]= float(MEMS[i].replace(',','.')) 
            time[i]= float(time[i].replace(',','.'))
        return np.array(Piezo), np.array(MEMS), np.array(time)
    else:
        for i in range(len(Piezo)):
            Piezo[i]= float(Piezo[i].replace(',','.'))   
            MEMS[i]= float(MEMS[i].replace(',','.')) 
            time[i]= float(time[i].replace(',','.'))
            Cap[i]= float(Cap[i].replace(',','.')) 
        return np.array(Piezo), np.array(MEMS), np.array(time), np.array(Cap)
        
def split_array(Piezo , MEMS, time, Cap):
    P, M, t, C = [], [], [], []
    POC_list, POC_ind = [], []
    X_val, Y_val = [], []
    X_ind = []
    
    for ind, elem in enumerate(Piezo):
        if elem == 'POC':
            POC_list.append(Piezo[ind+1])
            POC_ind.append(ind+1)
        if elem == 'X-Pos [um]':
            X_val.append(Piezo[ind+1])
            Y_val.append(MEMS[ind+1])
            X_ind.append(ind)
            
    for i,e in enumerate(POC_ind):
        if i==(len(POC_ind)-1):
            P.append(Piezo[POC_ind[i]+1:])        
            M.append(MEMS[POC_ind[i]+1:])
            t.append(time[POC_ind[i]+1:])
            C.append(Cap[POC_ind[i]+1:])
            
        else:
            P.append(Piezo[POC_ind[i]+1:X_ind[i+1]])        
            M.append(MEMS[POC_ind[i]+1:X_ind[i+1]])
            t.append(time[POC_ind[i]+1:X_ind[i+1]])
            C.append(Cap[POC_ind[i]+1:X_ind[i+1]])
        
    return P, M, t, C, POC_list, X_val, Y_val     

    
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

    for index_ul in range(index_h, len(Piezo_np)):
        piezo_del = Piezo_np[index_ul]-Piezo_np[index_ul-1]
        if abs(piezo_del) < 0.3:
            break
    
    if index_ul>np.argmin(MEMS_np):
        index_ul = np.argmin(MEMS_np)
        
    return [index_l, index_h, index_ul]
        
    
    

    # for index_ul in range(index_h,len(MEMS_np)):
    #     mems_del = MEMS_np[index_ul]-MEMS_np[index_ul-1]
    #     if mems_del > 2:
    #         break


 

def area_sphere(h_c, R=7500):
    A = np.pi*(2*R*h_c-h_c**2)
    return A

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
    return curve_fit(fit_func, reversed_piezo[start_index:end_index], reversed_MEMS[start_index:end_index], *p0, maxfev=100000,method='trf')

def detect_poc(F, h, delta=3):
    for index, value in enumerate(F):
        if np.abs(value) > 0+delta:
            return index
        
def poc_detect(Piezo_np, MEMS_np, time_np, Cap_np=['str'], delta = 2):
    index_poc = np.argmin(MEMS_np[0:np.argmax(MEMS_np)])
    if Cap_np[0]=='str':
        return Piezo_np[index_poc:], MEMS_np[index_poc:], time_np[index_poc:], index_poc
    else:
        return Piezo_np[index_poc:], MEMS_np[index_poc:], time_np[index_poc:],Cap_np[index_poc:], index_poc
        
    
    
    # for index_poc, val in enumerate(MEMS_np):
    #     if np.abs(val) > 0+delta:
    #         return Piezo_np[index_poc:], MEMS_np[index_poc:], time_np[index_poc:], Cap_np[index_poc:], index_poc
    #         break
        
def calc_stiff(curve_param, h_max):
    S = curve_param[1]*curve_param[0]*(h_max-curve_param[2])**(curve_param[1]-1)
    return S

def calc_hc(h_max, P_max, S, eps=0.762):
    h_s = eps*(P_max/S)
    h_c = h_max - h_s
    return h_c

def calc_hc_sph(Depth, Force, P_in):
    h_f = Depth[find_nearest(Force, P_in)]
    h_c = (np.max(Depth)+h_f)/2
    return h_c

def calc_emod(S, A, beta=1.05, nu_s=0.5, E_t=71, nu_t=0.17):
    # S in [nN/nm], A in [nm^2]-> E in [nN/nm^2]
    # daher: E[nN/nm^2]*10^9= E[Pa]
    E_r = 10**9 *(S*np.sqrt(np.pi))/(2*beta*(np.sqrt(A)))
    E = (1- nu_s**2)/(1/E_r-(1-nu_t**2)/(E_t*10**9))
    return E, E_r

def calc_H(P_max, A):
    return (P_max*10**-3)/(A*10**-12)

#Hertz Analysis

def calc_hertz(P, h, nu_s=0.5, E_t=71, nu_t=0.17):
    popt,pcov = curve_fit(func_hertz, h, P)
    E = (1- nu_s**2)/(1/popt-(1-nu_t**2)/(E_t*10**9))  
    return popt, E

def find_nearest(Force, value=0):
    delta = np.abs(Force-value)
    return np.argmin(delta)
    
def calc_JKRp(Depth, Force, R=7500, nu_s=0.5, E_t=71, nu_t=0.17):
    P_adh = np.min(Force)
    delta_adh = Depth[np.argmin(Force)]
    index_st = np.argmax(Force)
    index_end = np.argmin(Force)
    delta_0 = Depth[index_st + find_nearest(Force[index_st:index_end])]
    E_r = 10**9* ((-3*P_adh)/np.sqrt(R)) * ((3*(delta_0-delta_adh))/(1+4**(-2/3)))**-(3/2) 
    E = (1- nu_s**2)/(1/E_r-(1-nu_t**2)/(E_t*10**9))
    return E_r, E

def JKR_fitting_lmfit(params,P,h):
    a_R = params['a_R']
    h_contact = params['h_contact']
    P_adh = params['P_adh']
    h_fit = (a_R) * ((1+np.sqrt(1-P/P_adh)))**(4/3) -((2/3)*a_R) * ((1+np.sqrt(1-P/P_adh))/2)**(1/3) + h_contact
    return h_fit-h

def calc_a0(a_R, R=7500):
    return np.sqrt(a_R*R)
    
def calc_delta(P_adh, R=7500):
    return -(P_adh*2)/(3*np.pi*R)

def calc_jkr_Er(delta, a_0, R=7500):
    E_r = (9*np.pi*R**2*delta)/(2*a_0**3)
    return E_r*10**3

def JKR_analysis(P,h,params,R=7500, nu_s=0.5, E_t=71, nu_t=0.17):
    fitted_params = minimize(JKR_fitting_lmfit, params, args=(P, h), method='leastsq')
    a_R = fitted_params.params['a_R'].value
    h_c = fitted_params.params['h_contact'].value 
    P_adh = fitted_params.params['P_adh'].value
    a_0 = np.sqrt(a_R*R)
    delta = -(P_adh*2)/(3*np.pi*R)
    E_r = 10**3*(9*np.pi*R**2*delta)/(2*a_0**3)
    E = E = (1- nu_s**2)/(1/E_r-(1-nu_t**2)/(E_t*10**9))
    return [E, E_r, delta, a_0, P_adh, fitted_params]

def JKR_fit1(P, a_R, h_contact, P_adh):
    return (a_R) * ((1+np.sqrt(1-P/P_adh)))**(4/3) -((2/3)*a_R) * ((1+np.sqrt(1-P/P_adh))/2)**(1/3) + h_contact
  

#Für Kalibrierung de Kapazitäts-Verschiebung Koeffizienten
def calc_S(Cap_np_raw, Piezo_np_raw, index_l):
    S = (Cap_np_raw[index_l-1]-Cap_np_raw[0])/(Piezo_np_raw[index_l-1]-Piezo_np_raw[0])
    return S

def Hertz(P):
    E=6*10**(-9)
    R=7500
    return ((3*P)/(4*E*R**0.5))**(2/3)
    