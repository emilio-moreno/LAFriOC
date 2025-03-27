import os
import pandas as pd

date = '27-02-25'
dirpath = f'./PixelFly/{date}/'
listdir = os.listdir(dirpath)
saturated = []

def time_to_float(string):
	if 'us' in string:
		return float(string.replace('us', '')) / 1000
	return float(string.replace('ms', ''))


# Attenuation factor is actually 1E-4, but I'll
# divide by 1E-1 and convert from mW to uW (for the remaining 1E-3 attenuation)
attenuation_factor = 1E-1
print(attenuation_factor)
for file in listdir:
	if '.png' in file or 'X' in file: continue
	if "S" in file:
		power, exposition = file.split('_')[0], file.split('_')[1]
		power, exposition = float(power.replace('mW', '')) * attenuation_factor,  time_to_float(exposition.replace('S', ''))
		saturated.append((power, exposition))

df = pd.DataFrame(data=saturated, columns=["Power(uW)", "Exposure(ms)"])
print(df)
df.to_csv(os.path.join(dirpath, f"saturation_{date}.csv"))
