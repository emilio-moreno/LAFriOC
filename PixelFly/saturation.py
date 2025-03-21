import matplotlib.pyplot as plt
import pandas as pd
import sys, os
import numpy as np
from scipy.optimize import curve_fit
import json
sys.path.insert(0, '../../MiniPys/formatter')
import minipy_formatter as MF
MF.Format().rcUpdate()


def main_loop(data, ax, all_ax = None):
	df = pd.read_csv(data['saturation']['csvpath'])
	power, exposure = df['Power(uW)'], df['Exposure(ms)']

	all_power.extend(power.to_list())
	all_exposure.extend(exposure.to_list())


	figpath = curvefit_plot(power, exposure, data['date'], data['filter_factor'], ax, all_ax)

def curvefit_plot(power, exposure, date, filter_factor, ax, all_ax = None):
	# Curve fit
	if plot_inverse:
		f = lambda x, a: a * x	
		power = 1 / power
		xlabel = 'Power$^{-1}$ (uW$^{-1}$)'
		fitlabel = "Fit: $a x$"
	else:
		f = lambda x, a: a / x
		xlabel = 'Power (uW)'
		fitlabel = "Fit: $a / x$"
	opt, cov = curve_fit(f, power, exposure)
	opt, cov = float(opt[0]), float(cov[0, 0])

	# Plots
	P = np.linspace(min(power), max(power), 1000)
	for axis in [ax, all_ax]:
		if not axis: continue
		
		axis.scatter(power, exposure, label=f'Data - {filter_factor}%')
		fitlabel = f"{fitlabel}\n$a={opt:.2E}\\pm{cov:.2E}$"
		axis.plot(P, f(P, opt), label=fitlabel)

		axis.set(xlabel=xlabel, ylabel='Exposure (ms)')
		axis.legend(loc='upper right')
		axis.grid(True)
	ax.set_title(f"{date} / 780nm, Filter factor = {filter_factor}%")


def main():
	global all_power, all_exposure
	all_power, all_exposure = [], []

	global plot_inverse
	plot_inverse = False

	extra_str = ''
	if plot_inverse: extra_str = 'inverse_'

	# I/O
	with open('metadata.json', 'r') as f:
		global metadata
		metadata = json.load(f)

	mosaic = [[filter_name for filter_name in metadata]] + [len(metadata) * ['ALL']]
	fig, axs = plt.subplot_mosaic(mosaic, figsize=(18, 13), dpi=300)
	fig.suptitle("PixelFly - Maximum exposure before saturation")

	for filter_name in metadata:
		main_loop(metadata[filter_name], axs[filter_name], axs['ALL'])
	axs['ALL'].axhline(780, ls='--', color='red', label='Thermal saturation: 780ms')
	axs['ALL'].legend(loc='upper right')
	axs['ALL'].set_title("Summary")

	plt.subplots_adjust(hspace=0.33, wspace=0.33)
	figpath = f"./Saturation/{extra_str}saturation_780nm_27-02-25+13-03-25_E-2,4,6%.png"
	plt.savefig(figpath)
	#plt.show()
	
	#with open('metadata.json', 'w') as f:
	#	metadata = json.dump(metadata, f, indent=4)

if __name__ == '__main__':
	main()
