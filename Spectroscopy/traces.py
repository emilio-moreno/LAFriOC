import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import json


# %% Format
import os, sys
MF_formatter = r"C:\Users\dell 5810\Documents\M\GitHub\MiniPys\Formatter"
sys.path.insert(1, MF_formatter)
import minipy_formatter as MF
global Format
Format = MF.Format()
Format.rcUpdate()

os.chdir(r"C:\Users\dell 5810\Documents\M\Espectroscopía")


# %% Function definitions

def read_csv(filename):
    df = pd.read_csv(filename, usecols=[3, 4], names=['time', 'voltage'])
    return df


def smooth_array(x_values, array, avg_size):
    kernel = np.array(np.ones(avg_size)) / avg_size
    smoothed = np.convolve(kernel, array, mode='same')
    smoothed = smoothed[avg_size:-avg_size]
    x_values = x_values[avg_size:-avg_size]
    return x_values, smoothed
   

def format_plot(axs, filename, well):
    plt.suptitle(f"Pozo: {well} - {filename}")
    for ax in axs:
        axs[ax].set(title="Trazas", xlabel='Time', ylabel='Amplitude')
        axs[ax].grid(True)
        axs[ax].legend(loc='upper right')
    
    plt.subplots_adjust(hspace=0.75)
    
def get_peaks(time, voltage, prominence):
    if tmax != None and tmin != None:
        minfilter = time > tmin
        time, voltage = time[minfilter], voltage[minfilter]
        maxfilter = time < tmax
        time, voltage = time[maxfilter], voltage[maxfilter]
        
    peak_time = []
    peak_voltage = []
    found_peaks = find_peaks(voltage, prominence=prominence)[0]
    for index in found_peaks:
        peak_time.append(time[index])
        peak_voltage.append(voltage[index])
    
    return peak_time, peak_voltage


def mark_peaks(axs, peak_time, peak_voltage):
    axs['original'].scatter(peak_time, peak_voltage, s=7, color='red',
                            zorder=3)
    

def mainloop(filename):
    # I/O
    df = read_csv(filename)
    time, voltage = np.array(df['time']), np.array(df['voltage'])
    time, voltage = smooth_array(time, voltage, avg_size)
    # derivative = np.gradient(voltage)
    
    # Plotting
    # mosaic = [['original'], ['derivative']]
    mosaic = [['original']]
    fig, axs = plt.subplot_mosaic(mosaic, dpi=200, figsize=(7, 5))
    
    axs['original'].plot(time, voltage, label='Smoothed data')
    # axs['derivative'].plot(time, derivative, label='Derivative')
    peak_time, peak_voltage = get_peaks(time, voltage, prominence)
    mark_peaks(axs, peak_time, peak_voltage)   
    
    format_plot(axs, filename, well)     
    

def load_metadata(metadata_path):
    global metadata
    with open(metadata_path, 'r') as f:
        metadata = json.load(f)
    
    
def read_metadata(filename):
    global avg_size, prominence, well
    data = metadata[filename]
    avg_size, well = data['avg_size'], data['well']
    prominence= data['prominence']
    
    global tmax, tmin
    tmax, tmin = None, None
    if 'tmax' and 'tmin' in data:
        tmin, tmax = data['tmin'], data['tmax']
    

# %% I/O
dirpath = '18-03-25/20250311_trazas_espectrocopia'
load_metadata(os.path.join(dirpath, 'metadata_18-03-25.json'))
listdir = os.listdir(dirpath)

for filename in listdir:
    if not 'TEK' in filename: continue
    read_metadata(filename)
    mainloop(os.path.join(dirpath, filename))
