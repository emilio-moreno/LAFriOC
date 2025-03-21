import matplotlib.pyplot as plt
import pandas as pd
import sys, os
import numpy as np
from scipy.optimize import curve_fit
import json
sys.path.insert(0, '../../MiniPys/formatter')
import minipy_formatter as MF
MF.Format().rcUpdate()

colors = MF.Colors(('miku', '#599e94'), ('darkmiku', '#466964'), ('rin', '#ecd38c'),
                    ('darkrin', '#c6b176'), ('lgray', '#c0c0c0'))

#with open('metadata.json', 'r') as f:
#	metadata = json.load(f)

#data = metadata['0%']
dirpath = './Examples'
listdir = os.listdir(dirpath)
print(listdir)
#df = pd.read_csv(os.path.join(dirpath, listdir[1]))
#print(df.describe())

# maxfilt = df['Index'] < 16E3
# minfilt = df['Index'] > 15E3 
# index = df['Index'][maxfilt][minfilt]
# value = df[' Value'][maxfilt][minfilt]

# fig, ax = plt.subplots()

# ax.scatter(index, value, color=colors['miku'], s=5)
# #plt.fill(index, value, color=colors['miku'], alpha=0.2)
# #ax.scatter(index, value, color='purple')
# plt.show()
