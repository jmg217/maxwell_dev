#!/usr/bin/env_python

import matplotlib.pyplot as plt
import pandas as pd


data = pd.read_csv('results100.txt', sep="\t", header = None)

x= data['X.3']
time=data['X.13']


plt.ylabel('Runtime (sec)')
plt.xlabel('Mesh Parameter (b)')

plt.plot(x,time,'r', linewidth=2)

plt.savefig('timeMeshParam.pdf')
