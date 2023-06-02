import numpy as np
import matplotlib.pyplot as plt

truesigma1= np.array([228., 228., 227., 226., 226., 225., 224., 224., 223., 222., 222.,
       220., 220., 219., 218., 217., 216., 215., 213., 212., 211., 209.,
       131., 122., 117., 112., 107., 103.,  99.,  96.,  92.,  88.,  84.,
        80.,  77.,  72.,  68.,  63.,  58.,  52.,  44.,  33.])
truesigma2= np.array([150., 160., 167., 174., 183.])
truesigma3= np.array([207., 205., 203., 199., 193.])

temps1= np.array([86.        , 86.04081633, 86.08163265, 86.12244898, 86.16326531,
       86.20408163, 86.24489796, 86.28571429, 86.32653061, 86.36734694,
       86.40816327, 86.44897959, 86.48979592, 86.53061224, 86.57142857,
       86.6122449 , 86.65306122, 86.69387755, 86.73469388, 86.7755102 ,
       86.81632653, 86.85714286, 86.89795918, 86.93877551, 86.97959184,
       87.02040816, 87.06122449, 87.10204082, 87.14285714, 87.18367347,
       87.2244898 , 87.26530612, 87.30612245, 87.34693878, 87.3877551 ,
       87.42857143, 87.46938776, 87.51020408, 87.55102041, 87.59183673,
       87.63265306, 87.67346939])
temps2 =np.array([86.89795918, 86.93877551, 86.97959184, 87.02040816, 87.06122449])
temps2= temps2[::-1]
temps3 =np.array([86.89795918, 86.93877551, 86.97959184, 87.02040816, 87.06122449])

#find the index where the gradient of truesgima1 is most negative
splitIndex1=np.argmin(np.gradient(truesigma1))

#split truesigma1 and temps1 into two arrays at the index where the gradient is most negative
truesigma1a=truesigma1[:splitIndex1]
truesigma1b=truesigma1[splitIndex1:]
temps1a=temps1[:splitIndex1]
temps1b=temps1[splitIndex1:]

#join the arrays in the following order truesigma1a, truesigma3, truesigma2, truesigma1b
truesigma=np.concatenate((truesigma1a,truesigma3,truesigma2,truesigma1b))
temps=np.concatenate((temps1a,temps3,temps2,temps1b))

#order the values of truesigma by descending order
truesigma=truesigma[np.argsort(truesigma)]
#reverse the order of the values of truesigma
truesigma=truesigma[::-1]

plt.plot(temps,truesigma,'#1f77b4')
# plt.plot(temps1,truesigma1,'#1f77b4')
# plt.plot(temps2,truesigma2,'#1f77b4')
# plt.plot(temps3,truesigma3,'#1f77b4')

plt.show()

