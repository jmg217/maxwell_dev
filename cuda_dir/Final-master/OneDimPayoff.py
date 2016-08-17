#!/usr/bin/env_python

import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

data = pd.read_csv('OneDimHighPayoff.txt', sep="\t", header = None)

fig=plt.figure()
ax=fig.gca(projection='3d')

ax.set_xlabel('Time (Years)')
ax.set_ylabel('Stock Price ($)')
ax.set_zlabel('Option Value ($)')
plt.title('Evolution of High Estimator Option Values')

x=data['X.1']
z=data['X.2']
y=data['X.3']

ax.scatter(x,y,z,c=z)
#ax.plot_wireframe(x,y,z)
#fig.colorbar(surf,shrink=0.5, aspect=5)

plt.savefig('OptValEvol.pdf')
#ax.plot(x,y,z)


