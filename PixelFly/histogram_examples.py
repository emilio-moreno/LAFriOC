import matplotlib.pyplot as plt
import pandas as pd
import sys, os
import cv2 as cv
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
img = [file for file in listdir if ".png" in file and not "desc" in file]
# desc = [file for file in listdir if ".png" in file and "desc" in file]
txt = [file for file in listdir if ".txt" in file]

for txtpath, imgpath in zip(txt, img): #, desc):
	image = cv.imread(os.path.join(dirpath, imgpath))
	#description = cv.imread(os.path.join(dirpath, descpath))
	df = pd.read_csv(os.path.join(dirpath, txtpath))
	maxfilt1 = df['Index'] < 250
	minfilt1 = df['Index'] > 175
	maxfilt2 = df['Index'] < 16384
	minfilt2 = df['Index'] > 1000

	# index1 = df['Index'][minfilt1][maxfilt1]
	# value1 = df[' Value'][minfilt1][maxfilt1]
	index2 = df['Index'][minfilt2][maxfilt2]
	value2 = df[' Value'][minfilt2][maxfilt2]

	mosaic = [['hist2', 'img']]
	labels = ['Intensidades para saturación', #'Intensidades para saturación',
			  'Fotografía de la PixelFly'] #, 'Descripción']
	fig, axs = plt.subplot_mosaic(mosaic, figsize=(30, 10), dpi=300)
	title = txtpath.replace('.txt', '').replace('S', '').replace('_', ', ')
	fig.suptitle(title)
	
	# axs['hist1'].scatter(index1, value1, color=colors['miku'], s=5)
	axs['hist2'].scatter(index2, value2, color=colors['miku'], s=5)
	axs['hist2'].axvline(16384, color='red')
	axs['img'].imshow(image)
	#axs['desc'].imshow(description)

	for ax_label, label in zip(axs, labels):
		axs[ax_label].set_title(label)
		if "hist" in ax_label:
			axs[ax_label].set(xlabel='Bits', ylabel='Cuentas (#)')
			axs[ax_label].grid(True)

	plt.savefig("Examples/example" + imgpath)


# #plt.fill(index, value, color=colors['miku'], alpha=0.2)
# #ax.scatter(index, value, color='purple')
# plt.show()