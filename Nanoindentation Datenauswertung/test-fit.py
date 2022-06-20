# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 09:59:24 2022

@author: eim01
"""

import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

fit_range=0.5
alpha = 0.3
m = 1.2
h_f = 40
noise = np.random.normal(0,1,101)

h = np.linspace(0, 100, 101)
P = alpha*(h-h_f)**m
P_noise = alpha*(h-h_f)**m+noise

plt.subplot(2,1,1)
plt.plot(h, P, label='ideal curve')
plt.plot(h, P_noise, label='added noise')
plt.legend()


def func_exp(reversed_piezo, alpha, m, h_f):
    return alpha*(reversed_piezo - h_f)**m

def func_log(reversed_piezo, alpha, m, h_f):
    return np.log(alpha)+m*np.log(reversed_piezo-h_f)
    
def func_lin(reversed_piezo, alpha, m, h_f):
    return alpha+m*(reversed_piezo-h_f)

param_l, cov_l = curve_fit(func_lin, h[int(fit_range*len(h)):], P_noise[int(fit_range*len(h)):])
Y_l = func_lin(h, param_l[0], param_l[1], param_l[2])

param_exp, cov_exp = curve_fit(func_exp, h[int(fit_range*len(h)):], P_noise[int(fit_range*len(h)):], (0.5,1.3,1))
Y_exp = func_exp(h, param_exp[0], param_exp[1], param_exp[2])

param_log, cov_log = curve_fit(func_log, h[int(fit_range*len(h)):], np.log(P_noise[int(fit_range*len(h)):]), (0.5,1.3,1))
Y_log = np.exp(func_log(h, param_log[0], param_log[1], param_log[2]))



plt.subplot(2,1,2)
plt.plot(h, P_noise, label='data')
plt.plot(h, Y_l, label='fitted linear curve')
#plt.plot(h, Y_exp, label='power law fit')
plt.plot(h, Y_log, label='log fit')
plt.legend()