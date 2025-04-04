import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
import json
from hyperfine import HyperfineStructure


# %% Format
import os, sys
MF_formatter = r"C:\Users\dell 5810\Documents\M\GitHub\MiniPys\Formatter"
sys.path.insert(1, MF_formatter)
import minipy_formatter as MF
global Format
Format = MF.Format()
Format.rcUpdate()

os.chdir(r"C:\Users\dell 5810\Documents\M\GitHub\LAFriOC\Spectroscopy")


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
   

def format_plot(axs, filename, well, rescale=False):
    plt.suptitle(f"Pozo: {well} - {filename}")
    for ax in axs:
        if rescale: xlabel = "Frequency (MHz)"
        else: xlabel = "Time"
        axs[ax].set(title="Traces", xlabel=xlabel, ylabel="Amplitude")
        axs[ax].grid(True)
        axs[ax].legend(loc='upper right')
    
    plt.subplots_adjust(hspace=0.75)
    
def get_peaks(time, voltage, prominence, ignore_t=False):
    if tmax != None and tmin != None and not ignore_t:
        minfilter = time > tmin
        time, voltage = time[minfilter], voltage[minfilter]
        maxfilter = time < tmax
        time, voltage = time[maxfilter], voltage[maxfilter]
    
    global peak_time, peak_voltage
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
    

def mainloop(filename, rescale=False):
    # I/O
    df = read_csv(filename)
    time, voltage = np.array(df['time']), np.array(df['voltage'])
    time, voltage = smooth_array(time, voltage, avg_size)
    
    peak_time, peak_voltage = get_peaks(time, voltage, prominence)
    if rescale:
        scale = set_scale()
        time = time * scale
        peak_time, peak_voltage = get_peaks(time, voltage, prominence)
    # derivative = np.gradient(voltage)
    
    # Plotting
    # mosaic = [['original'], ['derivative']]
    mosaic = [['original']]
    fig, axs = plt.subplot_mosaic(mosaic, dpi=200, figsize=(7, 5))
    
    axs['original'].plot(time, voltage, label='Smoothed data')
    # axs['derivative'].plot(time, derivative, label='Derivative')
    mark_peaks(axs, peak_time, peak_voltage)   
    
    format_plot(axs, filename, well, rescale=True)     
    

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
        

def set_scale():
    # Using Rb87 P3/2 F3 as peak
    peak_diff_f = HyperfineStructure.third_crossover
    peak_diff_t = np.abs(peak_time[1] - peak_time[2])
    scale = peak_diff_f / peak_diff_t
    # peak_0_freq = peak_0_time * scale
    # zero = peak_0_freq - 193.7408
    return scale
    
        
def main():
    dirpath = '18-03-25/20250311_trazas_espectrocopia'
    load_metadata(os.path.join(dirpath, 'metadata_18-03-25.json'))
    listdir = os.listdir(dirpath)

    for filename in listdir:
        if not 'TEK' in filename: continue
        read_metadata(filename)
        mainloop(os.path.join(dirpath, filename))
        

def main_current():
    os.chdir(r"C:\Users\dell 5810\Documents\M\GitHub\LAFriOC\Spectroscopy\18-03-25")

    dirpath = r".\20250318_trazas_corriente"
    listdir = os.listdir(dirpath)

    load_metadata(os.path.join(dirpath, "metadata.json"))
    for directory in listdir:
        if not "ALL" in directory: continue
        
        for filename in os.listdir(os.path.join(dirpath, directory)):
            if not "CH1" in filename: continue
            if "CH1" in filename: continue
            read_metadata(filename)
            mainloop(os.path.join(dirpath, directory, filename),
                     rescale=True)
            

#%% majinai
if __name__ == '__main__':
    main_current()
    

