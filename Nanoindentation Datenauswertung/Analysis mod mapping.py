# -*- coding: utf-8 -*-
"""
Created on Thu Jun  29 13:00:06 2022

@author: eim01
"""

from functions import *
from scipy.optimize import curve_fit

poc = 0     #point of contact
area_func = []
e_tip = 0
nu_tip = 0

fit_range = [0, 1]  #


#%%
#############################################################################
#%%
#partial unload+ mod mapping
#data import and conversion
path = 'data/Test Measurement/array-partial unload'
#path = 'data/Saphir/PM-Sa-150ms, 3um, 20nms ,3cycles'
Piezo_, MEMS_, time_, Cap_ = imp_data(path)
Piezo, MEMS, time, Cap, POC, X_val, Y_val = split_array(Piezo_, MEMS_, time_, Cap_)
#print(POC)


Results = []
indent_depth = []
fit_param = []

for i, e in enumerate(Piezo):
    
    P, M, t = data_conversion(Piezo[i], MEMS[i], time[i])
    # Piezo offset position and conversion to [nm]
    P = (P - P[0])*1000
    # detect POC
    P, M, t, poc = poc_detect(P, M, t)  
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
    while index[-1] < np.argmin(M):
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
path = 'data/PDMS/AR-PDMS-60-10-60-1,2um-3x3-1'   
path = 'data/PDMS/AR-PDMS-60-10-60-1,2um-4x4-1'

path = 'data/PDMS/AR-PDMS-125-10-125-2,5um-diff_LF-7x7-1'
path = 'data/PDMS/AR-PDMS-125-10-125-2,5um-diff_LF-5x5-2'
path = ['O:/5-1/5-11/Messungen/2022/06_Nico_MA/05_Datenauswertung/Python/Nanoindentation Datenauswertung/data/PDMS/AR-SA-600-1-600-4um-diff_LF-5x5-1']
path = ['data/AR-SA-200-1-200-4um-diff_LF-4x4-1']
path = ['data/AR-SA-600-1-600-4um-diff_LF-4x4-1']
path = ['data/AR-SA-600-1-600-3,5um-diff_LF-4x4-2']
path = ['data/AR-SA-335-1-335-2,5um-same_LF-5x5-2']
path = ['data/AR-SA-300-1-300-2um-diff_LF-4x4-3']
#ab hier Messungen mit neuem MEMS
path = ['data/Saphir/AR-SA-450-1-450-3um-diff_LF-4x4-1']
path = ['data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-1']
path = ['data/Saphir/AR-SA-450-1-450-3um-diff_LF-4x4-2']
#path = ['data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-2']
#path = ['data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-3']
#path =['O:/5-1/5-11/Messungen/2022/06_Nico_MA/05_Datenauswertung/Python/Nanoindentation Datenauswertung/data/PDMS/AR-SA-600-1-600-4um-diff_LF-5x5-1','data/AR-SA-200-1-200-4um-diff_LF-4x4-1','data/AR-SA-600-1-600-4um-diff_LF-4x4-1', 'data/AR-SA-600-1-600-3,5um-diff_LF-4x4-2', 'data/AR-SA-300-1-300-2um-diff_LF-4x4-3']
#path = ['data/Saphir/AR-SA-450-1-450-3um-diff_LF-4x4-1', 'data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-1', 'data/Saphir/AR-SA-450-1-450-3um-diff_LF-4x4-2' ,'data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-2', 'data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-3']   #Messungen mit neuen MEMS
#path = ['data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-1','data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-2', 'data/Saphir/AR-SA-450-1-450-3um-diff_LF-5x5-3']   #Messungen mit neuen MEMS
bad_curves = 0
#S_m = 419*10**-6
P_offg, P_ing = [], []
S_l, S_ul, h_m = [], [], []
colours=['r','b','g']
for a, j in enumerate(path):
    Piezo_, MEMS_, time_, Cap_ = imp_data(j)
    Piezo, MEMS, time, Cap, POC, X_val, Y_val = split_array(Piezo_, MEMS_, time_, Cap_)
    S_load = []
    S_uload = []
    hmax = []
    P_in = []
    P_off = []
    n = 0
    x = 1
    s=0
    # if a==3:
    #     s=3
    # else:
    #     s=0

    for i in range(len(Piezo)):
               
        P, M, t, C = data_conversion(Piezo[i*x+n], MEMS[i*x+n], time[i*x+n], Cap[i*x+n])
        #M = C/S_m
        # Piezo offset position and conversion to [nm]
        # detect POC
        P_, M_, t_, C_, poc = poc_detect(P, M, t, C)  
        P = (P-P[0] )*1000
        P_c =  P-P[poc]
        # store each measurement as np array in a list
        Piezo[i], MEMS[i], time[i] = P, M, t
                
        Force = M * K
        Depth = (P - M)
        Depth = Depth - Depth[np.argmin(Force[0:np.argmax(Force)])]
        # plt.subplot(2,1,2)
        # plt.plot(M, Force, label = str(i))
        # plt.grid(b=True)
        # plt.xlabel('Eindringtiefe[nm]')
        # plt.ylabel('Kraft[nN]') 
        # plt.legend()
         
        index_l, index_h, index_ul = data_splitting(P, M, t)
        
        unload_Depth = Depth[index_h : index_ul+1]
        unload_Force = Force[index_h : index_ul+1] 
    
        reversed_Depth = unload_Depth[::-1] #[nm]
        reversed_Force = unload_Force[::-1] #[nN]
        
        unload_Piezo = P[index_h : index_ul+1]
        unload_MEMS = M[index_h : index_ul+1] 
    
        reversed_Piezo = unload_Piezo[::-1] #[nm]
        reversed_MEMS = unload_MEMS[::-1] #[nN]
        
        #plot segment boundaries
        index = [index_l, index_h, index_ul]
        #calculations
        
        #popt_exp, pcov_exp = fitting(reversed_Piezo, reversed_MEMS, fit_range, fit_func=func_exp)
        popt_load, pcov_load = curve_fit(func_lin, P[poc:index_l], M[poc:index_l])
        #popt_uload, pcov_uload = curve_fit(func_lin, reversed_Piezo, reversed_MEMS)
        
        popt_uload, pcov_uload = fitting(unload_Piezo[0:], unload_MEMS[0:],[0,1])
        
        if Force[index_ul]>-2600:
            bad_curves +=1
            S_uload.append(np.nan)
            S_load.append(np.nan)
            hmax.append(np.nan)
            P_in.append(np.nan)
            P_off.append(np.nan)
            continue
        P = P + popt_load[1]
        S_uload.append(popt_uload[0])
        S_load.append(popt_load[0])
        hmax.append(np.max(M))
        P_in.append(Force[poc])
        P_off.append(Force[index_ul])
        P_ing.append(Force[poc])
        P_offg.append(Force[index_ul])
        
        # plt.figure()
        # plt.plot(P, M)
        # plt.plot(np.take(P, index), np.take(M, index), ls = '', marker = "o")
        # plt.xlabel('Piezo[nm]')
        # plt.ylabel('MEMS[nm]')
        # plt.grid(visible=True)
        # plt.title(path)  
    
    S_l.append(S_load)
    S_ul.append(S_uload)
    h_m.append(hmax)
    
    P_in=np.array(P_in)
    P_off=np.array(P_off)
    hmax=np.array(hmax)
    S_load = np.array(S_load)
    S_uload = np.array(S_uload)
    hmax = hmax[~np.isnan(S_uload)]
    S_uload = S_uload[~np.isnan(S_uload)]
    S_load = S_load[~np.isnan(S_load)]
    P_in = P_in[~np.isnan(P_in)]
    P_off = P_off[~np.isnan(P_off)]
    print('Messreihe'+str(a+1)+'/n'+'mean +- std:')
    print(np.mean(S_uload),np.std(S_uload))
    
    plt.subplot(2,2,1)
    plt.subplot(2,1,2)
    #plt.subplot(1,1,1)
    #plt.plot(hmax, S_load, marker = '.',label='load segment'+str(a))
    plt.plot(hmax, S_uload, colours[a], marker = '.', label='Messreihe'+ str(a+1), alpha=0.5)
    plt.grid(visible=True)
    plt.legend()
    
    # plt.subplot(2,1,2)
    # plt.plot(hmax, P_off, marker = '.', label ='Messung'+str(a))
    # #plt.plot(hmax, P_in, marker = '.', label ='snap-in')
    # plt.grid(visible=True)
    # plt.xlabel('Z-max [nm]')
    # plt.ylabel('pull-off [nN]') 
    # plt.legend()
    
S_ul=np.array(S_ul)
h_m=np.array(h_m)
h_mean =[]
S_mean = []
S_std = []
for i in range(len(S_ul[0])):
    l = S_ul[:,i]
    k =h_m[:,i]
    n,m=[],[]
    for j in range(len(l)):
        if np.isnan(l[j]):
            continue
        n.append(l[j])
        m.append(k[j])
    h_mean.append(np.mean(m))
    S_mean.append(np.mean(n))
    S_std.append(np.std(n))

#plot stiffness versus indentation depth with errorbars 
#plt.subplot(1,1,1)       
plt.errorbar(h_mean, S_mean, yerr = S_std, color= 'k', capsize=3, marker='x', label='$ S_{mean} \pm \sigma $ ')
plt.xlabel('MEMS-Verschiebung [nm]')
plt.ylabel('Steigung')
#plt.title('S als Mittelwert mit 1 $\sigma$ Standarabweichung ')
plt.ylim(1.002,1.012)
plt.xlim(0.0,2400)
plt.legend(fontsize=9)

plt.subplot(2,2,1)
plt.plot(P,M, 'k')
plt.xlabel('Piezo-Verschiebung [nm]')
plt.ylabel('MEMS-Verschiebung [nm]')
plt.grid(visible=True)
plt.subplot(2,2,2)
plt.plot(t,M, 'k')
plt.xlabel('Zeit [s]')
plt.grid(visible=True)

# P_offg =np.array(P_offg)
# P_offg = P_offg[P_offg<-2000]
# P_ing =np.array(P_ing)
#P_ing = P_ing[P_ing<-2000]

fig, (ax1, ax2) = plt.subplots(1,2)
ax1.hist(P_offg, 25,histtype='stepfilled', color='r', stacked=None)
ax1.set_ylabel('Anzahl')
ax1.set_xlabel('Pull-off Kraft [nN]',fontsize=14)

ax2.hist(P_ing, 25, histtype='bar', color='g', stacked=None)
ax2.set_xlabel('Snap-in Kraft [nN]',fontsize=14)

fig.suptitle('Spitze-Probe Interaktionen (Saphir-Probe)', fontsize=16)
#%%

for i in range(len(Piezo)):
    plt.figure()
    idx = np.argmax(MEMS[i])
    plt.plot(Piezo[i][0:idx], MEMS[i][0:idx]*3.01,'r', marker='.',markersize=2, linestyle='-')
    plt.xlabel('Piezopositzion [nm]', size=18)
    plt.ylabel('Kraft [nN]', size=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # plt.ylim(-200,5)
    # plt.xlim(P[0],50)
    plt.grid()
#%%
#this script is for analysis of array Measurement on PDMS

path = ['data/PDMS/15.09/AR-PDMS-100ms,4um,20nms,1,5um offset,27x1-1']
path = ['data/PDMS_10-1/19.09/AR-PDMS10-1-100ms,4,5um,20nms,2,5um offset,16x1(x5)']
path = ['C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS_10-1/30.09/AR-PDMS10-1-100ms,4um,20nms,2um offset,16x1(x3)']
#PDMS 10:1 sample, with distance between array points
path = ['C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS_10-1/30.09/AR-PDMS10-1-100ms,4um,20nms,2um offset,16x1(x3)']
#path =['C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS_10-1/30.09/AR-PDMS10-1-100ms,4um,40nms,2um offset,16x1(x3)']
#path =['C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS_10-1/30.09/AR-PDMS10-1-100ms,4um,40nms,2um offset,16x1(x3)']

P_offg, P_ing = [], []

for a, j in enumerate(path):
    Piezo_, MEMS_, time_, Cap_ = imp_data(j)
    Piezo, MEMS, time, Cap, POC, X_val, Y_val = split_array(Piezo_, MEMS_, time_, Cap_)
    F, D, rD,rF = [0]*len(Piezo), [0]*len(Piezo),[0]*len(Piezo),[0]*len(Piezo)
    S_load, S_uload = [], []
    hmax = []
    P_in, P_off = [], []
    E_r_jkr, E_jkr = [], []
    E_Op = [[],[]]
    jkr_results = np.ones(len(Piezo))
    S_ges = []
    pow_params = []
    n = 0
    x = 1
    s = 0
    for i in range(len(Piezo)-s):
               
        P, M, t, C = data_conversion(Piezo[i*x+n], MEMS[i*x+n], time[i*x+n], Cap[i*x+n])
        P_, M_, t_, C_, poc = poc_detect(P, M, t, C)  
        P = (P )*1000
    
        # store each measurement as np array in a list
        Piezo[i], MEMS[i], time[i] = P, M, t
                
        Force = M * K
        Depth = (P - M)
        Depth = Depth - Depth[np.argmin(Force[0:np.argmax(Force)])]
        F[i] = Force
        D[i] = Depth
        plt.subplot(3,1,1)
        plt.plot(P, M, label = str(i))
        plt.grid(visible=True)
        plt.xlabel('Piezo[nm]')
        plt.ylabel('MEMS [nm]') 
        plt.legend()
         
        index_l, index_h, index_ul = data_splitting(P, M, t)
        
        unload_Depth = Depth[index_h : index_ul+1]
        unload_Force = Force[index_h : index_ul+1] 
    
        reversed_Depth = unload_Depth[::-1] #[nm]
        reversed_Force = unload_Force[::-1] #[nN]
        
        rD[i], rF[i] = reversed_Depth,reversed_Force
        
        #plot segment boundaries
        index = [index_l, index_h, index_ul]
        #calculations

        hmax.append(np.max(Depth))
        P_in.append(Force[poc])
        P_off.append(Force[index_ul])
        P_ing.append(Force[poc])
        P_offg.append(Force[index_ul])
        
        #JKR fitting
        # params = Parameters()

        # params.add('a_R',value=2000, min=0,vary=True)
        # params.add('h_contact', value=2000,  vary=True)
        # params.add('P_adh', value=np.min(Force), max=np.min(Force), vary=True)
        # jkr_results[i] = JKR_analysis(reversed_Force, reversed_Depth, params)[0]

        
        frange = [0.7, 0.95]
        start_index = int(len(reversed_Force)*frange[0])
        idx_0 = find_nearest(reversed_Force,value=3)
        reversed_Force= reversed_Force-Force[index_ul]
        
        if reversed_Force[start_index]<0:
            frange[0] = idx_0/len(reversed_Force)
            print(frange) 
        
        popt_exp, pcov_exp = fitting(reversed_Depth, reversed_Force, frange, fit_func=func_exp)
        reversed_Force= reversed_Force+Force[index_ul]
        S = calc_stiff(popt_exp, reversed_Depth[-1])    #[nN/nm]
        h_c = calc_hc_sph(reversed_Depth, reversed_Force, P_in[i]) #[nm]
        #h_c = calc_hc(reversed_Depth[-1], reversed_Force[-1], S, eps=0.774) #[nm]
        E_Op[0].append(10**-6*calc_emod(S, area_sphere(h_c))[0])
        E_Op[1].append(10**-6*calc_emod(S, area_sphere(h_c))[1])
        pow_params.append(popt_exp)
        S_ges.append(S)
        E_r_jkr.append(10**-6*calc_JKRp(Depth, Force, R= 7500)[0])        
        E_jkr.append(10**-6*calc_JKRp(Depth, Force, R= 7500)[1])   
        
        plt.subplot(3,1,2)
        plt.plot(Depth, Force)
        plt.plot(np.take(Depth, index), np.take(Force, index), ls = '', marker = "o")
        plt.grid(visible=True)
    
    

    i_num = 16
    i_rep = 5
    E_mean, E_std, h, S_mean, S_std = [], [], [],[], []
    E_Op_mean, E_Op_std = [], []
    for c in range(i_num):
        E_mean.append(np.mean(E_jkr[c::i_num]))
        E_std.append(np.std(E_jkr[c::i_num]))
        h.append(np.mean(hmax[c::i_num]))
        S_mean.append(np.mean(S_ges[c::i_num]))
        S_std.append(np.std(S_ges[c::i_num]))
        E_Op_mean.append(np.mean(E_Op[1][c::i_num]))
        E_Op_std.append(np.std(E_Op[1][c::i_num]))
        
    fig, (ax1,ax2) = plt.subplots(2,1)
    ax1.errorbar(h, E_mean, yerr = E_std,capsize=3, c = 'k')
    ax1.scatter(hmax, E_jkr, c= 'r')
    plt.grid(visible=True)
    ax1.set_xlabel('Eindringtiefe [nm]',fontsize=14)
    ax1.set_ylabel('E-Modul [MPa]', fontsize=14) 
    
    ax2.scatter(hmax, S_ges, c= 'r')
    ax2.errorbar(h, S_mean, yerr = S_std, capsize=3, c='r', label='Steifigkeit')
    ax2.set_ylabel('Steifigkeit [nN/nm]')
    ax2.legend()
    ax3 = ax2.twinx()
    ax3.scatter(hmax, E_Op[1])
    ax3.errorbar(h, E_Op_mean, yerr = E_Op_std, capsize=3, c = 'g', label ='E-OP')
    ax3.set_ylabel('E_Op[MPa]')
    ax3.legend()
    plt.title(path)

P_offg =np.array(P_offg)
P_ing =np.array(P_ing)

fig, (ax1, ax2) = plt.subplots(1,2)
ax1.hist(P_offg, 10,histtype='stepfilled', color='r', stacked=None)
ax1.set_ylabel('Anzahl')
ax1.set_xlabel('Pull-off Kraft [nN]',fontsize=14)

ax2.hist(P_ing, 10, histtype='bar', color='g', stacked=None)
ax2.set_xlabel('Snap-in Kraft [nN]',fontsize=14)
fig.suptitle('Spitze-Probe Interaktionen', fontsize=16)

#%%
#just JKR fitting and analysis
#PDMS reference sample
path = ['C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS/15.09/AR-PDMS-100ms,4um,20nms,1,5um offset,27x1-1']

#PDMS 10:1, repeated measurements at the same point
#path=['C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS_10-1/19.09/AR-PDMS10-1-100ms,4,5um,20nms,2,5um offset,16x1(x5)']

#PDMS 10:1 sample, with distance between array points
path = ['C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS_10-1/30.09/AR-PDMS10-1-100ms,4um,20nms,2um offset,16x1(x3)']
#path =['C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS_10-1/30.09/AR-PDMS10-1-100ms,4um,40nms,2um offset,16x1(x3)']



#path = ['C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS_10-1/30.09/AR-PDMS10-1-100ms,4um,20nms,2um offset,16x1(x3)','C:/Users/nicoe/Spyder Projekte/Nanoindentation Analysis/Python-Nanoindentation-Analysis/Nanoindentation Datenauswertung/data/PDMS_10-1/30.09/AR-PDMS10-1-100ms,4um,40nms,2um offset,16x1(x3)' ]
plt.rcParams.update({'font.size': 14})
Material =('A','B')
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
colors = ['b','r']
Laderate = ['20 nm/s','40 nm/s']
for a, j in enumerate(path):
    P_offg, P_ing = [], []
    hmax, S_ges = [],[]
    Piezo_, MEMS_, time_, Cap_ = imp_data(j)
    Piezo, MEMS, time, Cap, POC, X_val, Y_val = split_array(Piezo_, MEMS_, time_, Cap_)
    F, D, rD,rF = [0]*len(Piezo), [0]*len(Piezo),[0]*len(Piezo),[0]*len(Piezo)
    jkr_results = []
    for i in range(len(Piezo)):
               
        P, M, t, C = data_conversion(Piezo[i], MEMS[i], time[i], Cap[i])
        P_, M_, t_, C_, poc = poc_detect(P, M, t, C)  
        P = (P -P[0])*1000
    
        # store each measurement as np array in a list
        Piezo[i], MEMS[i], time[i] = P, M, t
                
        Force = M * K
        Depth = (P - M)
        Depth = Depth - Depth[np.argmin(Force[0:np.argmax(Force)])]
        F[i] = Force
        D[i] = Depth
        # plt.subplot(3,1,1)
        # plt.plot(P, M, label = str(i))
        # plt.grid(visible=True)
        # plt.xlabel('Piezo[nm]')
        # plt.ylabel('MEMS [nm]') 
        # plt.legend()
         
        index_l, index_h, index_ul = data_splitting(P, M, t)
        
        unload_Depth = Depth[index_h : index_ul+1]
        unload_Force = Force[index_h : index_ul+1] 
    
        reversed_Depth = unload_Depth[::-1] #[nm]
        reversed_Force = unload_Force[::-1] #[nN]
        hmax.append(np.max(Depth))
        P_offg.append(np.min(Force))
        P_ing.append(np.min(Force[0:np.argmax(Force)]))
        rD[i], rF[i] = reversed_Depth,reversed_Force
        
        #plot segment boundaries
        index = [index_l, index_h, index_ul]
        #calculations

        
        #JKR fitting
        params = Parameters()
        hc_val = Depth[np.argmin(Force)]
        params.add('a_R',value=2000, min=0,vary=True)
        #params.add('h_contact', value=hc_val,min=hc_val,max=hc_val+300,  vary=True)
        params.add('h_contact', value=2000,  vary=True)
        params.add('P_adh', value=np.min(unload_Force), max=np.min(unload_Force), vary=True)
        #fit of loading segment
        #jkr_results.append(JKR_analysis(Force[poc:index_l], Depth[poc:index_l], params))
        #fit of unloading segment
        jkr_results.append(JKR_analysis(reversed_Force[100:], reversed_Depth[100:], params))
        a_R = jkr_results[i][5].params['a_R'].value
        h_c = jkr_results[i][5].params['h_contact'].value 
        P_adh = jkr_results[i][5].params['P_adh'].value 
        #Depth = Depth - h_c
        
        # frange = [0.7, 0.99]
        # start_index = int(len(reversed_Force)*frange[0])
        # idx_0 = find_nearest(reversed_Force,value=3)
        # reversed_Force= reversed_Force-Force[index_ul]
        
        # # if reversed_Force[start_index]<0:
        # #     frange[0] = idx_0/len(reversed_Force)
        # #     print(frange) 
         
        # popt_exp, pcov_exp = fitting(reversed_Depth, reversed_Force, frange, fit_func=func_exp)
        # reversed_Force= reversed_Force+Force[index_ul]
        # S = calc_stiff(popt_exp, reversed_Depth[-1])    #[nN/nm]
        # S_ges.append(S)
        #uncomment for visual inspection of fitting
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(Depth,Force, color='b',label='Messdaten',marker='x')
    #loading segment
    plt.plot(JKR_fit1(Force[poc:index_ul],a_R,h_c, P_adh), Force[poc:index_ul], color='r',label = 'JKR fit')
    #unloading segment
    #plt.plot(JKR_fit1(reversed_Force,a_R,h_c, P_adh)-h_c, reversed_Force, color='r',label = 'JKR fit')
    plt.plot([JKR_fit1(P_adh,a_R,h_c, P_adh)], [P_adh], color='k', marker='.', markersize=25,label = 'snap-out', linestyle='')
    plt.plot(Depth[np.argmin(Force[0:np.argmax(Force)])], Force[np.argmin(Force[0:np.argmax(Force)])] , color='g',linestyle='', marker='.', markersize=25,label = 'snap-in')
    plt.xlabel('Eindringtiefe [nm]')
    plt.ylabel('Kraft [nN]')
    plt.legend(fontsize=10)
    plt.grid(visible=True)
    plt.subplot(2,1,2)
    plt.plot(t,Force, color='b',marker='x')
    plt.xlabel('Zeit [s]')
    plt.ylabel('Kraft [nN]')
    plt.grid(visible=True)
 
         
        #plt.subplot(3,1,2)
        # plt.figure()
        # plt.plot(Depth, Force)
        # plt.plot(JKR_fit1(reversed_Force,a_R,h_c, P_adh), reversed_Force,  label = 'JKR fit')
        # #plt.plot(np.take(Depth, index), np.take(Force, index), ls = '', marker = "o")
        # plt.grid(visible=True)

    hmax = np.array(hmax)    
    jkr_results=np.array(jkr_results)
    E_jkr = jkr_results[:,0]  
    P_Adh = jkr_results[:,4] 
    delta = jkr_results[:,2] 
    #S_ges = np.array(S_ges)
    # plt.figure()
    # plt.plot(E_jkr)
    
    
    #has to be manually adjusted before analysis
    i_num = 16
    i_rep = 3
    E_mean, E_std, h, S_mean, S_std = [], [], [],[], []
    P_Adh_m, P_Adh_s = [], []
    for c in range(i_num):
        E_mean.append(np.mean(E_jkr[c::i_num]))
        E_std.append(np.std(E_jkr[c::i_num]))
        h.append(np.mean(hmax[c::i_num]))
        P_Adh_m.append(np.mean(P_Adh[c::i_num]))
        P_Adh_s.append(np.std(P_Adh[c::i_num]))
        #S_mean.append(np.mean(S_ges[c::i_num]))
        #S_std.append(np.std(S_ges[c::i_num]))
        
    h = np.array(h)
    ax1.errorbar(h, E_mean, yerr = E_std,capsize=4, linewidth=2 , c = colors[a], label = Laderate[a])
    #ax1.scatter(hmax, E_jkr, c= 'r')
    #ax1.set_ylim(0,10)
    ax1.grid(visible=True)
    ax1.set_xlabel('Eindringtiefe [nm]')
    ax1.set_ylabel('E-Modul [MPa]') 
    ax1.legend()
    
    #ax3.errorbar(h, S_mean, yerr = S_std,capsize=4, linewidth=2 , c = 'g', label = 'Steifigkeit')
    ax4 = ax3.twinx()
    ax4.plot(h, area_sphere(h), label ='Flächenfunktion')
    ax3.legend()
    
    ax2.errorbar(h, P_Adh_m, yerr = P_Adh_s, capsize=4, linewidth=2 , c = colors[a], label = Laderate[a])
    #ax2.scatter(hmax, P_offg, color=colors[a], label = 'aus Messdaten')
    ax2.grid(visible=True)
    ax2.set_xlabel('Eindringtiefe [nm]')
    ax2.set_ylabel('Adhäsionskraft [nN]')
    ax2.legend()
