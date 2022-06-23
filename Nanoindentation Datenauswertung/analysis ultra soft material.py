# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 11:59:50 2022

@author: eim01
"""
import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


gamma = 0
R = 0
P_off = (3/2)*gamma*np.pi*R


#DMT
nu_1 = 0.3
nu_2 = 0.3

R_dmt = 0.00001
E_1 = 10**8
E_2 = 10**11
h_dmt = np.linspace(0.00001,101)
W = 1
K = (4/3)*(1/((1-nu_1**2)/E_1+(1-nu_2**2)/E_2))
P_dmt = (K*np.sqrt(h_dmt*R_dmt)**3)/R_dmt - 2*np.pi*R_dmt*W

plt.plot(h_dmt, P_dmt)