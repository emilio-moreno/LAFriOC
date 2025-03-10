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


for file in listdir:
	if '.png' in file or 'X' in file: continue
	if "S" in file:
		power, exposition = file.split('_')[0], file.split('_')[1]
		saturated.append((float(power.replace('mW', '')), time_to_float(exposition.replace('S', ''))))

df = pd.DataFrame(data=saturated, columns=["Power(mW)", "Exposure(ms)"])
print(df)
df.to_csv(os.path.join(dirpath, f"saturation_{date}.csv"))
