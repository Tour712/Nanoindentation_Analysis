# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 13:44:25 2022

@author: eim01
"""

from pre_processing_plain import *


path = 'data/MEMS Resolution Test/MEMS step test,0,5nm,10Hz-1'
Piezo, MEMS, time = imp_data(path)
Piezo_np_raw, MEMS_np_raw, time_np_raw = data_conversion(Piezo, MEMS, time)

Piezo_np_raw = (Piezo_np_raw - Piezo_np_raw[0])*1000


plt.plot(time_np_raw, Piezo_np_raw,ls = '-', markersize = 4, label='Piezo')
plt.plot(time_np_raw, MEMS_np_raw,ls = '-', markersize = 4, label='MEMS')
plt.xlabel('Zeit [s]')
plt.ylabel('Verschiebung von Piezo/MEMS [nm]')
plt.title('MEMS Stufentest')
plt.legend()