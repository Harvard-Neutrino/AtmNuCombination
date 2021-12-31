import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np

class Digitalizer:
	def __init__(self, image, scale):
		self.image = img.imread(image)
		self.scale = img.imread(scale)
		self.image_x_size, self.image_y_size = self.image.shape[:2]


	def set_palette(self, scale_max, sclae_min, bins):
		self.palette_x = np.logspace(scale_min, scale_max, 2 * n_bins + 1)[1::2]
		self.palette_clr = []

		dy = self.image_y_size / bins
		for iy in range(bins)[::-1]:
		    c = self.scale[int(np.round((iy + 0.5) * dy)), int(np.round(0.5 * self.img_x_size))]
		    self.palette_clr.append(c)

	def digitalize(self, n_bins_x, n_bins_y):
		