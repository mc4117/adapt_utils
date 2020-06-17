#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 13:33:38 2020

@author: mc4117
"""

import pandas as pd
import numpy as np
import pylab as plt

df_real = pd.read_csv('final_result_nx4.0.csv')

df_test = pd.read_csv('final_result_nx1.5.csv')

df_test1 = pd.read_csv('final_result_nx1.0.csv')

df_test2 = pd.read_csv('final_result_nx0.75.csv')

df_test3 = pd.read_csv('final_result_nx0.5.csv')

df_test4 = pd.read_csv('final_result_nx0.25.csv')

error_list = []
error_list.append(0.0)
error_list.append(sum([(df_test['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test1['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test2['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test4['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
print(error_list)
plt.plot([0.5, 2/3, 1, 4/3, 2, 4], error_list)

plt.loglog([2/3, 1, 4/3, 2, 4], error_list, '-o')
plt.ylabel('Error norm (m)')
plt.xlabel(r'$\Delta x$ (m)')
plt.show()


logx = np.log([2/3, 1, 4/3, 2, 4])
log_error = np.log(error_list)
np.polyfit(logx, log_error, 1)

# plot

"""
import pylab as plt

plt.plot(df_real['x'], df_real['bath'], label = 'nx = 2')
plt.plot(df_test1['x'], df_test1['bath'], label = 'nx = 1')
plt.plot(df_test2['x'], df_test2['bath'], label = 'nx = 0.5')
plt.plot(df_test3['x'], df_test3['bath'], label = 'nx = 0.25')
plt.legend()
"""
