#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:40:52 2023

@author: seanbartz
"""

import numpy as np
import matplotlib.pyplot as plt

cross_lambda_6_muValues=np.array([0,
50,
100,
150,
200,
250,
275,
285,
290])

cross_lambda_6_Temps=np.array([119.1,
118.76,
116.5,
113.359,
109,
103.775,
100.875,
99.68,
99.069])

FO_lambda_6_muValues=np.array([295,
300,
350,
400,
450,
500,
600,
700,
800,
850,
900,
1000,
1050,
1100,
1200,
1300,
1400,
1500,
1600,
1700,
1800])

FO_lambda_6_Temps=np.array([98.45076531,
97.82857143,
91.33,
84.46,
77.67,
70.46,
56.81,
44.79,
34.18,
29.85,
26.05,
19.1,
16.68,
14.185,
10.49,
7.75,
5.689,
4.144,
2.936,
1.734,
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
525,
530])

cross_lambda_5_Temps=np.array([94.72,
94.1,
92.5,
89.46,
85.73,
80.95,
75.76,
69.63,
63.5,
57.45,
51.6,
48.78,
48.21825])

FO_lambda_5_muValues=np.array([535,
550,
600,
650,
670,
750,
800,
900,
1000,
1050,
1100,
1200,
1300,
1400,
1500,
1600,
1700])

FO_lambda_5_Temps=np.array([47.65994898,
46,
40.19,
33.64,
32.19,
26.53,
22.3,
16.97,
12.31,
10.75,
9.23,
6.96,
4.914,
3.432,
2.439,
1.471,
0])
cross_lambda_45_muValues=np.array([0,
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
595])

cross_lambda_45_Temps=np.array([75.99,
75.56,
74.14,
72.31,
69.1,
65.91,
62.14,
57.72,
52.61,
47.35,
42.33,
37.5,
35.05,
33.265])

FO_lambda_45_muValues=np.array([600,
650,
700,
800,
850,
900,
1000,
1050,
1100,
1200,
1300,
1400,
1500])

FO_lambda_45_Temps=np.array([32.82755102,
28.681,
24.887,
18.43,
15.694,
12.53,
9.63,
8.21,
7.016,
4.062,
3.198,
1.457,
0])

#define the standard colors
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#plot the data
plt.plot(cross_lambda_6_muValues,cross_lambda_6_Temps,'--',color=colors[0])
plt.plot(FO_lambda_6_muValues,FO_lambda_6_Temps,label='$\lambda_1=6$',color=colors[0])
plt.plot(cross_lambda_5_muValues,cross_lambda_5_Temps,'--',color=colors[1])
plt.plot(FO_lambda_5_muValues,FO_lambda_5_Temps,label='$\lambda_1=5$',color=colors[1])
plt.plot(cross_lambda_45_muValues,cross_lambda_45_Temps,'--',color=colors[2])
plt.plot(FO_lambda_45_muValues,FO_lambda_45_Temps,label='$\lambda_1=4.5$',color=colors[2])
plt.title('$m_q = 15$ MeV')

#make the plot look nice
plt.xlabel('$\mu$ (MeV)')
plt.ylabel('T (MeV)')
plt.xlim(0,1900)
plt.ylim(0,125)


#make a legend
plt.legend()

#save the plot as png, pdf, jpg, and eps
plt.savefig('phaseDiagramPlot_mq_15.png',dpi=500)
plt.savefig('phaseDiagramPlot_mq_15.pdf')
plt.savefig('phaseDiagramPlot_mq_15.jpg')
plt.savefig('phaseDiagramPlot_mq_15.eps')



plt.show()