import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import json
from traces import *


# %% Format
import os, sys
MF_formatter = r"C:\Users\dell 5810\Documents\M\GitHub\MiniPys\Formatter"
sys.path.insert(1, MF_formatter)
import minipy_formatter as MF
global Format
Format = MF.Format()
Format.rcUpdate()

os.chdir(r"C:\Users\dell 5810\Documents\M\GitHub\LAFriOC\Spectroscopy\18-03-25")

dirpath = r".\20250318_trazas_corriente"
listdir = os.listdir(dirpath)

metadata = load_metadata(os.path.join(dirpath, "metadata.json"))
for directory in listdir:
    if not "ALL" in directory: continue
    
    for filename in os.listdir(os.path.join(dirpath, directory)):
        if not "CH1" in filename: continue
        print(filename)
    

