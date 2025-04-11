'''Analysis of how total counts increase with exposure time.'''
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from extract_saturation import time_to_float

#%%
dirpath = r"./13-03-25/E-4"
listdir = os.listdir(dirpath)
files = [filename for filename in listdir if '.txt' in filename]
files = [filename for filename in files if 'X_' not in filename]

for filename in files:
    # exposure = filename.split('_')[1].replace('S', '')
    # exposure = time_to_float(exposure) # ms
    df = pd.read_csv(os.path.join(dirpath, filename))
    total_counts = sum(df[' Value'])
    print(total_counts)