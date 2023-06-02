#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt



mu =np.array([39.5,
	185,
	302,
	530,
	685,
	895,
	1195])


T = np.array([78.106,
	77.842,
	70.07065,
	48.21825,
	34.67079592,
	21.3455,
	10.85])

mq=np.array([3,6,9,15,19,24,30])

# normalize y values between 0 and 1
mq_normalized = (mq - np.min(mq)) / (np.max(mq) - np.min(mq))

# define colormap
viridis = plt.cm.get_cmap('viridis')

# Create the figure and axis objects
fig, ax = plt.subplots()

# plot the line
for i in range(len(mu)-1):
    ax.scatter(mu[i:i+2], T[i:i+2],color=viridis(mq_normalized[i]))

# Add the legend
sm = plt.cm.ScalarMappable(cmap=viridis, norm=plt.Normalize(vmin=np.min(mq), vmax=np.max(mq)))
sm._A = []  # Set the dummy variable
cbar = plt.colorbar(sm)
cbar.ax.set_title('Quark Mass')

plt.ylabel('T (MeV)')
plt.xlabel('$\mu$ (MeV)')
plt.title('Critical points for $\lambda_1 = 5$')
plt.tight_layout()

plt.savefig('criticalPoints_lambda5.png',dpi=500)

# Show the plot
plt.show()


