import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

plt.rcParams.update({"font.family": "Times New Roman", "mathtext.fontset": "cm"})

filename = "wavelength_vs_temperature.csv"
df = pd.read_csv(filename)
df = df[~df['temperature(C)'].str.contains('\*')]
temperature = np.fromiter((t.replace('*', '') for t in df['temperature(C)']), dtype=np.float64)
frequency = np.fromiter((f.split('/')[0] for f in df['frequency(THz)']), dtype=np.float64)

frequency = frequency[temperature > 19.11]
temperature = temperature[temperature > 19.11]

fig, ax = plt.subplots()
ax.scatter(temperature, frequency, s=15, c='r')
ax.grid(True, linestyle='--', alpha=0.7)
ax.set(title='Frequency vs Temperature', xlabel='Temperature (C)', ylabel='Frequency (THz)')
plt.show()