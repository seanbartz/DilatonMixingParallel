#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:40:52 2023

@author: seanbartz
"""

import numpy as np
import matplotlib.pyplot as plt

cross_lambda_7438_muValues=np.array([0,
50,
100,
150,
200,
250,
300,
350,
400,
425,
432,
433])

cross_lambda_7438_Temps=np.array([151.11,
150.26,
148.89,
145.95,
141.97,
136.21,
130.58,
123.95,
117.11,
113.55,
112.5220838,
112.3726])

FO_lambda_7438_muValues=np.array([434,
435,
436,
450,
500,
550,
575,
580,
585,
595,
600])

FO_lambda_7438_Temps=np.array([112.22511,
112.075,
111.929,
109.74,
102.16,
94.47,
91.41633163,
90.57059559,
89.55435235,
87.75510204,
86.73])

cross_lambda_5_muValues=np.array([775,
800,
825,
835,
845,
850,
860,
865,
870,
875,
885,
890,
895])

cross_lambda_5_Temps=np.array([30,
27.5,
25.8,
25.1,
24.425,
24.32,
23.475,
23.15,
22.83,
22.5375,
21.9325,
21.6375,
21.3455])

cross_lambda_78_muValues=np.array([0,
25,
50,
75,
100,
150,
200,
205,
210])

cross_lambda_78_Temps=np.array([155.91,
155.76,
155.33,
154.61,
153.61,
150.785,
146.9225,
146.4828,
146.48275])

FO_lambda_78_muValues=np.array([215,
225,
250])

FO_lambda_78_Temps=np.array([145.57,
144.632,
141.12])

#plot the data
plt.plot(cross_lambda_7438_muValues,cross_lambda_7438_Temps,'--',label='$\lambda=7.438$')
plt.plot(FO_lambda_7438_muValues,FO_lambda_7438_Temps,label='$\lambda=7.438$')
plt.plot(cross_lambda_5_muValues,cross_lambda_5_Temps,'--',label='$\lambda=5$')
plt.plot(cross_lambda_78_muValues,cross_lambda_78_Temps,'--',label='$\lambda=7.8$')
plt.plot(FO_lambda_78_muValues,FO_lambda_78_Temps,label='$\lambda=7.8$')

#make a legend
plt.legend()
plt.show()