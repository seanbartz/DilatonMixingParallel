#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:40:52 2023

@author: seanbartz
"""

import numpy as np
import matplotlib.pyplot as plt

cross_lambda_24_muValues=np.array([0,
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

cross_lambda_24_Temps=np.array([100.36,
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

FO_lambda_24_muValues=np.array([900,
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

FO_lambda_24_Temps=np.array([21.0577551020408,
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

cross_lambda_15_muValues=np.array([0,
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

cross_lambda_15_Temps=np.array([94.72,
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

FO_lambda_15_muValues=np.array([535,
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

FO_lambda_15_Temps=np.array([47.65994898,
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

cross_lambda_9_muValues=np.array([0,
50,
100,
150,
200,
250,
300,
302])

cross_lambda_9_Temps=np.array([88.76,
88.49,
86.67,
83.84,
80.21,
75.49,
70.28653,
70.07065])

FO_lambda_9_muValues=np.array([303,
304,
305,
310,
325,
350,
400,
450,
500,
550,
600,
650,
700,
750,
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
1600])

FO_lambda_9_Temps=np.array([69.96229861,
69.85351,
69.744327,
69.1935,
67.5126,
64.6352,
58.6996,
52.92,
47.14,
41.243,
36.15,
31.23,
26.665,
22.753,
19.38,
16.379,
13.43,
10.04,
8.25,
7.05,
5.224,
3.865,
1.989,
0.851,
0])

plt.title('$\lambda=5$')
#define the standard colors
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.plot(cross_lambda_24_muValues,cross_lambda_24_Temps,'--',color=colors[0])
plt.plot(FO_lambda_24_muValues,FO_lambda_24_Temps,color=colors[0],label='$m_q=24$')
plt.plot(cross_lambda_15_muValues,cross_lambda_15_Temps,'--',color=colors[1])
plt.plot(FO_lambda_15_muValues,FO_lambda_15_Temps,label='$m_q=15$',color=colors[1])
plt.plot(cross_lambda_9_muValues,cross_lambda_9_Temps,'--',color=colors[2])
plt.plot(FO_lambda_9_muValues,FO_lambda_9_Temps,label='$m_q=9$',color=colors[2])
#label the axes
plt.xlabel('$\mu$ (MeV)')
plt.ylabel('$T$ (MeV)')

#make a legend
plt.legend()

#save plot as a pdf in a subfolder called "plots"
plt.savefig('plots/phaseDiagramPlot_mq_9_15_24_lambda_5.pdf')

#save plot as an eps in a subfolder called "plots"
plt.savefig('plots/phaseDiagramPlot_mq_9_15_24_lambda_5.eps')

#save plot as a jpg in a subfolder called "plots"
plt.savefig('plots/phaseDiagramPlot_mq_9_15_24_lambda_5.jpg')

plt.show()