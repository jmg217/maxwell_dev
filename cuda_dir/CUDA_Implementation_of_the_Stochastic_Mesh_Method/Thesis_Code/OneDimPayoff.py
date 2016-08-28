#!/usr/bin/env_python

import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab

data = pd.read_csv('OneDimHighPayoff.txt', sep="\t", header = None)

fig=plt.figure()
ax=fig.gca(projection='3d')

ax.set_ylabel('Time (Years)')
ax.set_xlabel('Stock Price ($)')
ax.set_zlabel('Option Value ($)')
plt.title('Evolution of High Estimator Put Option Values')

x=data['X.1']
z=data['X.2']
y=data['X.3']


ax.scatter(y,x,z,c=z)

#pylab.xlim([0,200])
#ax.plot_wireframe(x,y,z)
#fig.colorbar(surf,shrink=0.5, aspect=5)

#NOTE THIS CHANGES THE POSITION OF Z AXIS
tmp_planes = ax.zaxis._PLANES 
ax.zaxis._PLANES = ( tmp_planes[2], tmp_planes[3], 
                     tmp_planes[0], tmp_planes[1], 
                     tmp_planes[4], tmp_planes[5])
view_1 = (25, -135)
view_2 = (25, -45)
init_view = view_2
ax.view_init(elev=20, azim=320)
pylab.xlim([0,200])
plt.savefig('OptValEvolPut1.png')

#ax.plot(x,y,z)


