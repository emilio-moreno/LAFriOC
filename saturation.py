import matplotlib.pyplot as plt
import pandas as pd
import sys, os
import numpy as np
from scipy.optimize import curve_fit
sys.path.insert(0, '../MiniPys/formatter')
import minipy_formatter as MF
MF.Format().rcUpdate()

# I/O
date = '27-02-25'
dirpath = f"./PixelFly/{date}/"
filename = os.path.join(dirpath, f"saturation_{date}.csv")
df = pd.read_csv(filename)
power, exposure = df['Power(mW)'], df['Exposure(ms)']

# Curve fit
f = lambda x, a: a / x
opt, cov = curve_fit(f, power, exposure)

# Plots
P = np.linspace(min(power), max(power), 1000)
fig, axs = plt.subplot_mosaic([['Saturated']])
ax = axs['Saturated']
ax.scatter(power, exposure, label='Data')
ax.plot(P, f(P, opt[0]), label=f'Fit: $a / x$\n$a={float(opt[0]):.2f}\\pm{float(cov[0]):.5f}$')
ax.set(title=f'PixelFly ({date})\nMaximum exposure before saturation (780 nm)', xlabel='Power (mW)',
					 ylabel='Exposure (ms)')
ax.grid(True)
ax.legend()
plt.savefig(os.path.join(dirpath, f'saturation_780nm_{date}.png'))
plt.show()