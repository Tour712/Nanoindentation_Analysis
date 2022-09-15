# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 09:41:29 2022

@author: eim01
"""
import numpy as np
from lmfit import Parameters,minimize, fit_report

x = np.linspace(0,10,100)

# Creating random data for y axis 
# Here the slope (m) is predefined to be 2.39645
# The intercept is 0
# random.normal method just adds some noise to data
y = 2.39645 * x + np.random.normal(0, 2, 100)



# Define the fitting function
def linear_fitting_lmfit(params,x,y):
    m = params['m']
    c = params['c']
    y_fit = m*x + c
    return y_fit-y

# Defining the various parameters
params = Parameters()
# Slope is bounded between min value of 1.0 and max value of 3.0
params.add('m', min=1.0, max=3.0)
# Intercept is made fixed at 0.0 value
params.add('c', value=0.0, vary = False)

# Calling the minimize function. Args contains the x and y data.
fitted_params = minimize(linear_fitting_lmfit, params, args=(x,y,), method='least_squares')

# Getting the fitted values
m = fitted_params.params['m'].value
c = fitted_params.params['c'].value    

# Printing the fitted values
print('The slope (m) is ', m)
print('The intercept (c) is ', c)

# Pretty printing all the statistical data
print(fit_report(fitted_params))