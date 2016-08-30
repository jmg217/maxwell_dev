#!/usr/bin/env_python

import matplotlib.pyplot as plt
import pandas as pd


data = pd.read_csv('resultsDim.txt', sep="\t", header = None)

x= data['X.1']
time=data['X.14']


plt.ylabel('Time (sec)')
plt.xlabel('Number of Underlying Assets')
plt.title('Computational Cost of Multiple Underlying Assets')
plt.plot(x,time,'r', linewidth=2)

plt.savefig('Dim.pdf')

