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
433,])

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
600,
675,
775,
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
895,
900,
1000,
1250,
1400,
1500,
1650,
1800,
1900])

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
86.73,
75.83892617,
62.14285714,
58.7755102,
55.71428571,
54.48979592,
53.46938776,
52.55714286,
51.53673469,
50.63,
50.19,
49.47,
48.49,
47.87,
47.26,
46.85,
36.88,
19.96,
13.8,
10.67,
6.83,
3.23,
0])

cross_lambda_5_muValues=np.array([0,
50,
100,
150,
200,
250,
300,
350,
400,
450,
500,
550,
575,
580,
585,
595,
600,
675,
775,
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

cross_lambda_5_Temps=np.array([100.36,
99.14,
97.75,
94.75,
91.04,
86.5,
81.36,
75.52,
69.53,
63.52,
57.53,
51.47,
48.52,
47.9,
47.2,
46.06122449,
45.96326531,
39.40816327,
30,
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

FO_lambda_5_muValues=np.array([900,
1000,
1050,
1075,
1100,
1200,
1250,
1300,
1400,
1500,
1600,
1700,
1750])

FO_lambda_5_Temps=np.array([21.0577551020408,
16.2,
13.94,
13.08,
12.25,
9.23,
8.36,
7.08,
5.26,
3.99,
2.96,
1.23,
0])
cross_lambda_78_muValues=np.array([0,
25,
50,
75,
100,
150,
200,
205,
210,])

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
250,
300,
350,
400,
450,
500,
550,
600,
700,
800,
900,
1000,
1100,
1200,
1300,
1400,
1500,
1600,
1700,
1800,
1900,
2000,
2100])

FO_lambda_78_Temps=np.array([145.57,
144.632,
141.12,
136.435,
130.121,
123.26,
115.9,
108.31,
100.68,
92.76,
77.74,
63.39,
51.07,
40.3,
31.62,
24.94,
19.55,
15.11,
11.62,
9.23,
7.35,
5.21,
3.1,
1.9,
0])
#define the standard colors
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#plot the data
plt.plot(cross_lambda_7438_muValues,cross_lambda_7438_Temps,'--',color=colors[1])
plt.plot(FO_lambda_7438_muValues,FO_lambda_7438_Temps,label='$\lambda=7.438$',color=colors[1])
plt.plot(cross_lambda_5_muValues,cross_lambda_5_Temps,'--',color=colors[0])
plt.plot(FO_lambda_5_muValues,FO_lambda_5_Temps,label='$\lambda=5$',color=colors[0])
plt.plot(cross_lambda_78_muValues,cross_lambda_78_Temps,'--',color=colors[2])
plt.plot(FO_lambda_78_muValues,FO_lambda_78_Temps,label='$\lambda=7.8$',color=colors[2])

#make the plot look nice
plt.xlim(0,2100)
plt.ylim(0,160)
plt.xticks(np.arange(0,2001,200))
plt.yticks(np.arange(0,161,20))


#label the axes
plt.xlabel('$\mu$ (MeV)')
plt.ylabel('T (MeV)')

#title the plot with the value of the quark mass
plt.title('$m_q=24$ MeV')

#make a legend
plt.legend()

#save the plot as a pdf, jpg, png, and eps in a folder called "plots"
plt.savefig('plots/phaseDiagram_3LambdaVals_mq_24.pdf')
plt.savefig('plots/phaseDiagram_3LambdaVals_mq_24.jpg')
plt.savefig('plots/phaseDiagram_3LambdaVals_mq_24.png', dpi=500)
plt.savefig('plots/phaseDiagram_3LambdaVals_mq_24.eps')



plt.show()