#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 17:37:12 2020

@author: mc4117
"""

import pandas as pd
import pylab as plt
import numpy as np

df_real = pd.read_csv('fixed_output/bed_trench_output_uni_4.csv')
df_real1 = pd.read_csv('fixed_output/bed_trench_output_uni_2.0.csv')
df_real2 = pd.read_csv('fixed_output/bed_trench_output_uni_0.8.csv')
df_real3 = pd.read_csv('fixed_output/bed_trench_output_uni_0.4.csv')
df_real4 = pd.read_csv('fixed_output/bed_trench_output_uni_0.2.csv')
df_real4a = pd.read_csv('fixed_output/bed_trench_output_uni_0.125.csv')
df_real4b = pd.read_csv('fixed_output/bed_trench_output_uni_0.1.csv')
df_real5a = pd.read_csv('fixed_output/bed_trench_output_uni_0.05.csv')

error_list = []
#error_list.append(0.0)
error_list.append(sum([(df_real1['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_real2['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_real3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_real4['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_real4a['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_real4b['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_real5a['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
print(error_list)

plt.loglog([0.1, 0.25, 0.5, 1, 2], error_list, '-o')
plt.ylabel('Error norm (m)')
plt.xlabel(r'$\Delta x$ (m)')
plt.show()


logx = np.log([0.1, 0.25, 0.5, 1, 2])
log_error = np.log(error_list)
poly = np.polyfit(logx, log_error, 1)
print(poly[0])


data = pd.read_csv('experimental_data.csv', header=None)
df_exp0 = pd.read_csv('fixed_output/bed_trench_output4.csv')
df_exp = pd.read_csv('fixed_output/bed_trench_output2.0.csv')
df_exp1 = pd.read_csv('fixed_output/bed_trench_output0.8.csv')
df_exp2 = pd.read_csv('fixed_output/bed_trench_output0.4.csv')
df_exp3 = pd.read_csv('fixed_output/bed_trench_output0.2.csv')
df_exp4 = pd.read_csv('fixed_output/bed_trench_output0.1.csv')

error_listexp = []
error_listexp.append(np.sqrt(sum([(df_exp0['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp))])))
error_listexp.append(np.sqrt(sum([(df_exp['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp))])))
error_listexp.append(np.sqrt(sum([(df_exp1['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp))])))
error_listexp.append(np.sqrt(sum([(df_exp2['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp))])))
error_listexp.append(np.sqrt(sum([(df_exp3['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp))])))
error_listexp.append(np.sqrt(sum([(df_exp4['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp))])))

plt.loglog([0.05, 0.1, 0.25, 0.5, 1, 2], error_listexp, '-o')
plt.ylabel('Error norm (m)')
plt.xlabel(r'$\Delta x$ (m)')
plt.show()
