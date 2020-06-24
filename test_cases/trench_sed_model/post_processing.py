#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 17:37:12 2020

@author: mc4117
"""

import pandas as pd

df_real = pd.read_csv('fixed_output/bed_trench_output1.0.csv')
df_real1 = pd.read_csv('fixed_output/bed_trench_output0.8.csv')
df_real2 = pd.read_csv('fixed_output/bed_trench_output0.6.csv')
df_real3 = pd.read_csv('fixed_output/bed_trench_output0.4.csv')
df_real4 = pd.read_csv('fixed_output/bed_trench_output0.2.csv')

stop

df_real = pd.read_csv('fixed_output/bed_trench_output2.0.csv')
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
df_test3_3 = pd.read_csv('adapt_output2/bed_trench_output_0.5_25.0.csv')
df_test4_3 = pd.read_csv('adapt_output2/bed_trench_output_0.5_50.0.csv')
df_test5_3 = pd.read_csv('adapt_output2/bed_trench_output_0.5_100.0.csv')
df_test6_3 = pd.read_csv('adapt_output2/bed_trench_output_0.5_200.0.csv')
#df_test7_3 = pd.read_csv('adapt_output2/bed_trench_output_0.5_500.0.csv')

error_list3 = []
error_list3.append(sum([(df_test1_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list3.append(sum([(df_test2_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list3.append(sum([(df_test3_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list3.append(sum([(df_test4_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list3.append(sum([(df_test5_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list3.append(sum([(df_test6_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
#error_list3.append(sum([(df_test7_3['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))

df_test1_4 = pd.read_csv('adapt_output2/bed_trench_output_0.6_0.0.csv')
df_test2_4 = pd.read_csv('adapt_output2/bed_trench_output_0.6_10.0.csv')
df_test3_4 = pd.read_csv('adapt_output2/bed_trench_output_0.6_25.0.csv')
df_test4_4 = pd.read_csv('adapt_output2/bed_trench_output_0.6_50.0.csv')
df_test5_4 = pd.read_csv('adapt_output2/bed_trench_output_0.6_100.0.csv')
df_test6_4 = pd.read_csv('adapt_output2/bed_trench_output_0.6_200.0.csv')
#df_test7_4 = pd.read_csv('adapt_output2/bed_trench_output_0.6_500.0.csv')

error_list4 = []
error_list4.append(sum([(df_test1_4['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list4.append(sum([(df_test2_4['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list4.append(sum([(df_test3_4['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list4.append(sum([(df_test4_4['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list4.append(sum([(df_test5_4['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list4.append(sum([(df_test6_4['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
#error_list4.append(sum([(df_test7_4['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))


df_test1_5 = pd.read_csv('adapt_output2/bed_trench_output_0.8_0.0.csv')
df_test2_5 = pd.read_csv('adapt_output2/bed_trench_output_0.8_10.0.csv')
df_test3_5 = pd.read_csv('adapt_output2/bed_trench_output_0.8_25.0.csv')
df_test4_5 = pd.read_csv('adapt_output2/bed_trench_output_0.8_50.0.csv')
df_test5_5 = pd.read_csv('adapt_output2/bed_trench_output_0.8_100.0.csv')
df_test6_5 = pd.read_csv('adapt_output2/bed_trench_output_0.8_200.0.csv')
#df_test7_5 = pd.read_csv('adapt_output2/bed_trench_output_0.8_500.0.csv')

error_list5 = []
error_list5.append(sum([(df_test1_5['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list5.append(sum([(df_test2_5['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list5.append(sum([(df_test3_5['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list5.append(sum([(df_test4_5['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list5.append(sum([(df_test5_5['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list5.append(sum([(df_test6_5['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
#error_list5.append(sum([(df_test7_5['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))df_test1_4 = pd.read_csv('adapt_output2/bed_trench_output_0.6_0.0.csv')

df_test1_6 = pd.read_csv('adapt_output2/bed_trench_output_1.0_0.0.csv')
df_test2_6 = pd.read_csv('adapt_output2/bed_trench_output_1.0_10.0.csv')
df_test3_6 = pd.read_csv('adapt_output2/bed_trench_output_1.0_25.0.csv')
df_test4_6 = pd.read_csv('adapt_output2/bed_trench_output_1.0_50.0.csv')
df_test5_6 = pd.read_csv('adapt_output2/bed_trench_output_1.0_100.0.csv')
df_test6_6 = pd.read_csv('adapt_output2/bed_trench_output_1.0_200.0.csv')
#df_test7_5 = pd.read_csv('adapt_output2/bed_trench_output_0.8_500.0.csv')

error_list6 = []
error_list6.append(sum([(df_test1_6['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list6.append(sum([(df_test2_6['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list6.append(sum([(df_test3_6['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list6.append(sum([(df_test4_6['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list6.append(sum([(df_test5_6['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
error_list6.append(sum([(df_test6_6['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))
#error_list5.append(sum([(df_test7_5['bath'][i] - df_real['bath'][i])**2 for i in range(len(df_real))]))df_test1_4 = pd.read_csv('adapt_output2/bed_trench_output_0.6_0.0.csv')

