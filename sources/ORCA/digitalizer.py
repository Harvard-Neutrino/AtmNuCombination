import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np
from util import gaus_fit
from params import *

class Digitalizer:
	def __init__(self, image, scale):
		self.image = img.imread(image)
		self.scale = img.imread(scale)
		self.scale_y_size, self.scale_x_size = self.scale.shape[:2]
		self.image_y_size, self.image_x_size = self.image.shape[:2]
		self.extracted = None
		self.gaussians = None


	def set_palette(self, scale_max, scale_min, bins):
		self.palette_x = np.logspace(scale_min, scale_max, 2 * bins + 1)[1::2]
		self.palette_clr = []

		dy = self.scale_y_size / bins
		for iy in range(bins)[::-1]:
			c = self.scale[int(np.round((iy + 0.5) * dy)), int(np.round(0.5 * self.scale_x_size))]
			# print("in the ", iy, "th bin, the color extracted is ", c)
			self.palette_clr.append(c)

	def digitalize(self, n_bins_x, n_bins_y):
		ext_data = np.zeros((n_bins_x, n_bins_y))

		dx = self.image_x_size / n_bins_x
		dy = self.image_y_size / n_bins_y

		for ix in range(n_bins_x):       # loops over all pixels in the extracted data
		    for iy in range(n_bins_y):
		        c = self.image[int(np.round((iy + 0.5) * dy)), int(np.round((ix + 0.5) * dx))] # find color at center of this pixel
		        if abs(np.sum(c ** 2) - 4) < 1e-4: # Blank = 0
		            ext_data[ix, n_bins_y - 1 - iy] = 0
		        else:
		            ic = np.argmin(np.sum((c[None,:] - self.palette_clr) ** 2, axis = 1) )         # find closest color in palette
		            ext_data[ix, n_bins_y - 1 - iy] = self.palette_x[ic]

		self.extracted = ext_data

	def fit(self):
		self.gaussians = np.ndarray((self.extracted.shape[0],2), float)
		for i in range(self.extracted.shape[0]):
			sigma, mu = gaus_fit(self.extracted[i], x_bins)
			self.gaussians[i][0] = sigma
			self.gaussians[i][1] = mu