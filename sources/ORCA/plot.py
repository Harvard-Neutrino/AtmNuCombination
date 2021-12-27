import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np

## Digitize the scale ##
#scale_min = 0 # Minimum value of the scale
#scale_max = 0.01 * 183.42/25.5 # Maximum value of the scale
scale_min = -3 # Minimum value of the scale
scale_max = 0 # Maximum value of the scale
n_bins = 20 # Number of color bins

#palette_x = np.linspace(scale_min, scale_max, 2*n_bins+1)[1::2]
palette_x = np.logspace(scale_min, scale_max, 2*n_bins+1)[1::2]
palette_clr = []

img_data = img.imread("scale.png")
img_size_y, img_size_x = img_data.shape[:2]

dy = img_size_y / n_bins
for iy in range(n_bins)[::-1]:
    c = img_data[int(np.round((iy+0.5)*dy)), int(np.round(0.5*img_size_x))] # find color at center of this pixel
    palette_clr.append(c)

## Digitize the covariance matrix ##
n_bins_x  = 22
n_bins_y = 22

my_data = np.zeros((n_bins_x, n_bins_y))

#img_data = img.imread('migration.png')
img_data = img.imread('migrationCascades.png')
#img_data = img.imread('migrationTracks.png')
img_size_y, img_size_x = img_data.shape[:2]
dx = img_size_x / n_bins_x
dy = img_size_y / n_bins_y

for ix in range(n_bins_x):       # loops over all pixels in the extracted data
    for iy in range(n_bins_y):
        c = img_data[int(np.round((iy+0.5)*dy)), int(np.round((ix+0.5)*dx))] # find color at center of this pixel
        if abs(np.sum(c**2) - 4) < 1e-4: # Blank = 0
            my_data[ix, n_bins_y-1-iy] = 0
        else:
            ic = np.argmin( np.sum((c[None,:] - palette_clr)**2, axis=1) )         # find closest color in palette
            my_data[ix, n_bins_y-1-iy] = palette_x[ic]

#np.savetxt("Migration_matrixTracks.dat", my_data)
np.savetxt("Migration_matrixCascades.dat", my_data)

exit(0)
my_data = np.loadtxt("Migration_matrix.dat")
# If I multiply my_data by my true energy binned spectrum, I will get my reconstructed energy binned spectrum
print(np.matmul(my_data, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
print(np.matmul(my_data, [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
print(np.matmul(my_data, [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
print(np.matmul(my_data, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]))

epsilon = np.sum(my_data, axis=0)
d = np.loadtxt("Measured.dat")
u = np.loadtxt("MC_truth.dat")