#!/usr/bin/env_python

import matplotlib.pyplot as plt
import pandas as pd


data = pd.read_csv('results110five.txt', sep="\t", header = None)

x= data['X.2']
V=data['X.7']
v=data['X.5']
HighBound=data['X.10']
LowBound=data['X.9']
PointEst=data['X.11']

plt.ylabel('Option Price ($)')
plt.xlabel('Mesh Parameter (b)')
plt.title('Convergence of Mesh and Path Estimators')
plt.plot(x,V,'r', linewidth=2, label='Mesh Estimator')
plt.plot(x,v,'b', linewidth=2, label='Path Estimator')
plt.plot(x,LowBound,'b--', linewidth=2, label='Lower 90% Confidence Bound')
plt.plot(x,HighBound,'r--',linewidth=2, label='Upper 90% Confindence Bound')
plt.plot(x,PointEst,'g--',linewidth=2, label='Point Estimate')
plt.legend()
plt.savefig('figFive110two.pdf')

