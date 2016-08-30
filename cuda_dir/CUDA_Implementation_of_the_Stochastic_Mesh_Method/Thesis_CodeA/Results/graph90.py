#!/usr/bin/env_python

import matplotlib.pyplot as plt
#from plt import rc
import pandas as pd

#rc('text', usetex=True)

data = pd.read_csv('results90five.txt', sep="\t", header = None)



x= data['X.2']
V=data['X.7']
v=data['X.5']
HighBound=data['X.10']
LowBound=data['X.9']
PointEst=data['X.11']

plt.ylim([1,2])
plt.ylabel('Option Price ($)')
plt.xlabel('Mesh Parameter (b)')
plt.title('Convergence of Mesh and Path Estimators')
plt.plot(x,V,'r', linewidth=2, label='Mesh Estimator')
plt.plot(x,v,'b', linewidth=2, label='Path Estimator')
plt.plot(x,LowBound,'b--', linewidth=2, label='Lower 90% Confidence Bound')
plt.plot(x,HighBound,'r--',linewidth=2, label='Upper 90% Confindence Bound')
plt.plot(x,PointEst,'g--',linewidth=2, label='Point Estimate')
plt.legend()
plt.savefig('figFive90two.pdf')
