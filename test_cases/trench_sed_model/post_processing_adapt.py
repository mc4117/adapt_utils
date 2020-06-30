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
df_test1 = pd.read_csv('adapt_output/bed_trench_output_uni_0.1_1.csv')
df_test2 = pd.read_csv('adapt_output/bed_trench_output_uni_0.1_2.csv')
df_test3 = pd.read_csv('adapt_output/bed_trench_output_uni_0.1_3.csv')
df_test4 = pd.read_csv('adapt_output/bed_trench_output_uni_0.1_4.csv')
df_test5 = pd.read_csv('adapt_output/bed_trench_output_uni_0.1_5.csv')
df_test6 = pd.read_csv('adapt_output/bed_trench_output_uni_0.1_10.csv')
df_test7 = pd.read_csv('adapt_output/bed_trench_output_uni_0.1_12.csv')
df_test8 = pd.read_csv('adapt_output/bed_trench_output_uni_0.1_15.csv')
df_test9 = pd.read_csv('adapt_output/bed_trench_output_uni_0.1_17.csv')
df_test10 = pd.read_csv('adapt_output/bed_trench_output_uni_0.1_20.csv')
df_test11 = pd.read_csv('adapt_output/bed_trench_output_uni_0.1_40.csv')
df_test12 = pd.read_csv('adapt_output/bed_trench_output_uni_0.1_60.csv')


error_list = []
#error_list.append(0.0)
error_list.append(sum([(df_test1['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test2['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test4['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test5['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test6['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test7['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test8['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test9['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test10['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test11['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test12['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
print(error_list)

data = pd.read_csv('experimental_data.csv', header=None)

df_exp1 = pd.read_csv('adapt_output/bed_trench_output_0.1_1.csv')
df_exp2 = pd.read_csv('adapt_output/bed_trench_output_0.1_2.csv')
df_exp3 = pd.read_csv('adapt_output/bed_trench_output_0.1_3.csv')
df_exp4 = pd.read_csv('adapt_output/bed_trench_output_0.1_4.csv')
df_exp5 = pd.read_csv('adapt_output/bed_trench_output_0.1_5.csv')
df_exp6 = pd.read_csv('adapt_output/bed_trench_output_0.1_10.csv')
df_exp7 = pd.read_csv('adapt_output/bed_trench_output_0.1_12.csv')
df_exp8 = pd.read_csv('adapt_output/bed_trench_output_0.1_15.csv')
df_exp9 = pd.read_csv('adapt_output/bed_trench_output_0.1_17.csv')
df_exp10 = pd.read_csv('adapt_output/bed_trench_output_0.1_20.csv')
df_exp11 = pd.read_csv('adapt_output/bed_trench_output_0.1_40.csv')
df_exp12 = pd.read_csv('adapt_output/bed_trench_output_0.1_60.csv')

error_listexp = []
error_listexp.append(np.sqrt(sum([(df_exp1['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp1))])))
error_listexp.append(np.sqrt(sum([(df_exp2['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp1))])))
error_listexp.append(np.sqrt(sum([(df_exp3['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp1))])))
error_listexp.append(np.sqrt(sum([(df_exp4['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp1))])))
error_listexp.append(np.sqrt(sum([(df_exp5['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp1))])))
error_listexp.append(np.sqrt(sum([(df_exp6['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp1))])))
error_listexp.append(np.sqrt(sum([(df_exp7['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp1))])))
error_listexp.append(np.sqrt(sum([(df_exp8['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp1))])))
error_listexp.append(np.sqrt(sum([(df_exp9['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp1))])))
error_listexp.append(np.sqrt(sum([(df_exp10['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp1))])))
error_listexp.append(np.sqrt(sum([(df_exp11['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp1))])))
error_listexp.append(np.sqrt(sum([(df_exp12['bath'][i] - data[1].dropna()[i])**2 for i in range(len(df_exp1))])))

print(error_listexp)