import math
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import sterileprob
import sterileutil as su
import xarray as xr

def circlePlot():
	z = np.linspace(1/1000, np.pi  - 1/1000, 1000)
	x1= np.zeros(len(z))
	y1= np.zeros(len(z))
	x2= np.zeros(len(z))
	y2= np.zeros(len(z))
	for i in range(1, len(z)):
		x1[i] = sterileprob.prob(2, 1, "nova", "IH", "nova", m4 = 1, t14 = 0, t24= 0, t34= 0, \
						mand13 = True, vald13 = z[i], d14 = 0, d34 = 0, anti = False)
		y1[i] = sterileprob.prob(2, 1, "nova", "IH", "nova", m4 = 1, t14 = 0, t24= 0, t34= 0, \
						mand13 = True, vald13 = z[i], d14 = 0, d34 = 0, anti = True)
		x2[i] = sterileprob.prob(2, 1, "nova", "NH", "nova", m4 = 1, t14 = 0, t24= 0, t34= 0, \
						mand13 = True, vald13 = z[i], d14 = 0, d34 = 0, anti = False)
		y2[i] = sterileprob.prob(2, 1, "nova", "NH", "nova", m4 = 1, t14 = 0, t24= 0, t34= 0, \
						mand13 = True, vald13 = z[i], d14 = 0, d34 = 0, anti = True)
	plt.plot(x1,y1, color = "orange")
	plt.plot(x2,y2, color = "blue")
	plt.savefig("try.png")

def paperPlot():
	x = np.linspace(-1, 1, 1000)
	y = np.linspace(0.3, 0.7, 1000)
	x, y = np.meshgrid(x, y)
	z = np.zeros(len(x))
	z = sterileprob.prob(2, 1, 't2k', 'NH', 'nova', m4 = 1, t14 = np.arcsin(np.sqrt(y)), t24 = 0.15,\
		 t34= 0, mand13 = True, vald13 = x* np.pi,manE = False, valE = 0.9, d14 = 0, d34 = 0, anti = False, avg = False)

	t2kCI90 = 0.0675 + 1.645 * 0.0085

	plt.contour(x, y, z)
	plt.savefig("try2.png")

paperPlot()
# circlePlot()
