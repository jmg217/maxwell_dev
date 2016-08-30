#!/usr/bin/env_python

import matplotlib.pyplot as plt
import pandas as pd


data = pd.read_csv('stocksimulation.txt', sep="\t", header = None)

x= data.pop('X.1')

plt.ylabel('Stock Price ($)')
plt.xlabel('Time (Years)')
plt.title('Stock Price Simulation')
plt.plot(x,data)
plt.savefig('simulation.pdf')

