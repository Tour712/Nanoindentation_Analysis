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
            return Piezo_np[index_poc:], MEMS_np[index_poc:], time_np[index_poc:]
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
#%%
#partial unload
#data import and conversion
path = 'data/Test Measurement/array-partial unload'
Piezo_, MEMS_, time_, Cap_ = imp_data(path)
Piezo, MEMS, time, Cap, POC, X_val, Y_val = split_array(Piezo_, MEMS_, time_, Cap_)
print(POC)




# detect POC

Results = []
indent_depth = []
fit_param = []

for i, e in enumerate(Piezo):
    
    P, M, t = data_conversion(Piezo[i], MEMS[i], time[i])
    # Piezo offset position and conversion to [nm]
    P = (P - P[0])*1000
    # detect POC
    P, M, t = poc_detect(P, M, t)  
    # store each measurement as np array in a list
    Piezo[i], MEMS[i], time[i] = P, M, t
    
    #plt.plot(time[i], MEMS[i])
    
    #detect indices of load, hold and unload segment
    index, index_l, index_h, index_ul = [],[],[],[]
    index = data_splitting(P, M, t)
    index_l.append(data_splitting(P, M, t)[0])
    index_h.append(data_splitting(P, M, t)[1])
    index_ul.append(data_splitting(P, M, t)[2])
    
    #save end of segment index in list
    while index[-1] < (len(P)-1):
        ind = data_splitting(P, M, t, index_start= index[-1]+1)
        index = index + ind
        index_l.append(ind[0])
        index_h.append(ind[1])
        index_ul.append(ind[2])

    # fit curve and calculate Results
    unload_Piezo, unload_MEMS = [],[]
    popt_log, pcov_log = [], []
    S = []
    hmax = []
    
    for j in range(len(index_l)):
        up = P[index_h[j] : index_ul[j]]
        uM = M[index_h[j] : index_ul[j]]
        unload_Piezo.append(up[::-1])
        unload_MEMS.append(uM[::-1]) 
        hmax.append(up[0])
    
        
        if j==len(index_l):
            par, cov = fitting(unload_Piezo[j] , np.log(unload_MEMS[j]), [0.2, 0.95], (1.0, 1,0), fit_func=func_log)
            #uncomment for power law fit
            #par, cov = fitting(unload_Piezo[j] , np.log(unload_MEMS[j]), [0.3, 0.95], (1.0, 1.0,1), fit_func=func_exp)
        else:         
            par, cov = fitting(unload_Piezo[j] , np.log(unload_MEMS[j]), [0.3, 0.95], (1.0 ,1,0), fit_func=func_log) 
            #uncomment for power law fit
            #par, cov = fitting(unload_Piezo[j] , unload_MEMS[j], [0.3, 0.95], (1.0, 1.0,1), fit_func=func_exp)
        popt_log.append(par)
        pcov_log.append(cov)
        
       
        S.append(calc_stiff(par, up[0]))
    fit_param.append(popt_log) 
    indent_depth.append(hmax)
    Results.append(S)   #Results is a list of length(number of array points), where each entry is a list of lenght(number of load cycles)

#plot results
indent_depth = np.array(indent_depth)
Results = np.array(Results)
S_mean, S_std = [], []

for i in range(5):
    S_std.append(np.std(Results[:,i]))
    S_mean.append(np.mean(Results[:,i]))
    
#plot stiffness versus indentation depth with errorbars    
plt.errorbar(indent_depth[0], S_mean, yerr = S_std)
plt.xlabel('h [nm]')
plt.ylabel('Steifigkeit')


#     ax.plot(unload_Piezo[i], func_exp(unload_Piezo[i], par[0], par[1], par[2]), label = 'log fit' + str(i))


#%%

path = 'data/S_calib/S_calib_2308-2x2'
Piezo_, MEMS_, time_, Cap_ = imp_data(path)
Piezo, MEMS, time, Cap, POC, X_val, Y_val = split_array(Piezo_, MEMS_, time_, Cap_)
print(POC)

Results = []
indent_depth = []
fit_param = []

for i, e in enumerate(Piezo):
    
    P, M, t = data_conversion(Piezo[i], MEMS[i], time[i])
    # Piezo offset position and conversion to [nm]
    P = (P - P[0])*1000
    plt.plot(P, M)
    # detect POC
    P, M, t = poc_detect(P, M, t)  
    # store each measurement as np array in a list
    Piezo[i], MEMS[i], time[i] = P, M, t
    










