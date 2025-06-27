import matplotlib.pyplot as plt
import pandas as pd
import os

plt.rcParams.update({"font.family": "Times New Roman", "mathtext.fontset": "cm"})

filename = "power_vs_current_19-11C.csv"
df = pd.read_csv(filename)

fig, ax = plt.subplots()
ax.scatter(df['current(mA)'], df['power(mW)'], s=5, c='r')
ax.set(title='Power vs Current', xlabel='Current (mA)', ylabel='Power (mW)')
plt.show()