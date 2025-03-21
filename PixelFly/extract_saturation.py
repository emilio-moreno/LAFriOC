import os
import pandas as pd
import json

def time_to_float(string):
	if 'us' in string:
		return float(string.replace('us', '')) / 1000
	return float(string.replace('ms', ''))

def power_to_float(string, filter_factor):
	# As numbers are quite small I'll use uW
	power = float(string.replace('mW', '')) * 1E3 # uW
	# Filter factor is in percentage so I divide by 100.
	return power * filter_factor / 100

def main_loop(data):
	dirpath = data['dirpath']
	listdir = os.listdir(dirpath)
	filter_factor = 10**int(data['filter_factor'].replace('E', ''))
	saturated = []
	for file in listdir:
		if not "mW" in file or not '.txt' in file: continue
		# Excluded files.
		if "X" in file: continue
		# Using S on files to denote saturation.
		if "S" in file:
			power, exposition = file.split('_')[0], file.split('_')[1].replace('S', '')
			saturated.append((power_to_float(power, filter_factor), time_to_float(exposition)))

	df = pd.DataFrame(data=saturated, columns=["Power(uW)", "Exposure(ms)"])
	csvpath = f"./Saturation/saturation_{data['date']}_{data['filter_factor']}.csv"
	df.to_csv(csvpath)
	data["saturation"] = {"csvpath": csvpath}

	print(dirpath)
	print(df, "\n")

	return data


def main():
	# I/O
	with open('metadata.json', 'r') as f:
		metadata = json.load(f)
	for filename in metadata:
		metadata[filename] = main_loop(metadata[filename])

	with open('metadata.json', 'w') as f:
		metadata = json.dump(metadata, f, indent=4)


if __name__ == '__main__':
	main()

