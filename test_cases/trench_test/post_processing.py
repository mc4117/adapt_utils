#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 17:37:12 2020

@author: mc4117
"""

import pandas as pd

df_real = pd.read_csv('fixed_output/160_bed_trench_output.csv')
df_test1 = pd.read_csv('adapt_output2/bed_trench_output_0.2_0.0.csv')
df_test2 = pd.read_csv('adapt_output2/bed_trench_output_0.2_10.0.csv')
df_test3 = pd.read_csv('adapt_output2/bed_trench_output_0.2_25.0.csv')
df_test4 = pd.read_csv('adapt_output2/bed_trench_output_0.2_50.0.csv')
df_test5 = pd.read_csv('adapt_output2/bed_trench_output_0.2_100.0.csv')
df_test6 = pd.read_csv('adapt_output2/bed_trench_output_0.2_200.0.csv')
df_test7 = pd.read_csv('adapt_output2/bed_trench_output_0.2_500.0.csv')

error_list = []
error_list.append(sum([(df_test1['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test2['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test4['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test5['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test6['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list.append(sum([(df_test7['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))

df_test1_2 = pd.read_csv('adapt_output2/bed_trench_output_0.4_0.0.csv')
df_test2_2 = pd.read_csv('adapt_output2/bed_trench_output_0.4_10.0.csv')
df_test3_2 = pd.read_csv('adapt_output2/bed_trench_output_0.4_25.0.csv')
df_test4_2 = pd.read_csv('adapt_output2/bed_trench_output_0.4_50.0.csv')
df_test5_2 = pd.read_csv('adapt_output2/bed_trench_output_0.4_100.0.csv')
df_test6_2 = pd.read_csv('adapt_output2/bed_trench_output_0.4_200.0.csv')
df_test7_2 = pd.read_csv('adapt_output2/bed_trench_output_0.4_500.0.csv')

error_list2 = []
error_list2.append(sum([(df_test1_2['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list2.append(sum([(df_test2_2['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list2.append(sum([(df_test3_2['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list2.append(sum([(df_test4_2['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list2.append(sum([(df_test5_2['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list2.append(sum([(df_test6_2['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list2.append(sum([(df_test7_2['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))

df_test1_3 = pd.read_csv('adapt_output2/bed_trench_output_0.5_0.0.csv')
df_test2_3 = pd.read_csv('adapt_output2/bed_trench_output_0.5_10.0.csv')
#df_test3_3 = pd.read_csv('adapt_output2/bed_trench_output_0.5_25.0.csv')
df_test4_3 = pd.read_csv('adapt_output2/bed_trench_output_0.5_50.0.csv')
df_test5_3 = pd.read_csv('adapt_output2/bed_trench_output_0.5_100.0.csv')
df_test6_3 = pd.read_csv('adapt_output2/bed_trench_output_0.5_200.0.csv')
#df_test7_3 = pd.read_csv('adapt_output2/bed_trench_output_0.5_500.0.csv')

error_list3 = []
error_list3.append(sum([(df_test1_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list3.append(sum([(df_test2_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
#error_list3.append(sum([(df_test3_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list3.append(sum([(df_test4_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list3.append(sum([(df_test5_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list3.append(sum([(df_test6_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
#error_list3.append(sum([(df_test7_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))