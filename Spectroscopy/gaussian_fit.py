import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import json
from traces import read_csv, smooth_array

filename = r"C:\Users\dell 5810\Documents\M\GitHub\LAFriOC\Spectroscopy\18-03-25\20250318_trazas_corriente\ALL0004\F0004CH1.CSV"
df = read_csv(filename)

time, voltage = df['time'], df['voltage']
print(len(time), len(voltage))
time, voltage = smooth_array(time, voltage, avg_size=5)
plt.plot(time, voltage)


#%% Noise
f = lambda x, a, b: a * x + b
time_filter = time > 0.025
bg_time = time[time_filter]
bg_voltage = voltage[time_filter]

popt, cov = curve_fit(f, bg_time, bg_voltage)
plt.plot(time, voltage - f(time, *popt))

# plt.plot(time, f(time, *popt))