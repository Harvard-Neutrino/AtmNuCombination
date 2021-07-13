import h5py
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('fname', type=str, nargs='?', default='NULL')
args = parser.parse_args()
infile = args.fname

if infile == 'NULL':
    sys.exit('Please introduce your HDF5 file')

f1 = h5py.File(infile, 'r')
print(f1.keys())

for key in f1.keys():
#    print(key)
    ds = f1[key]
    data = np.array(ds[()])
#    print(data)
    plt.hist(data, bins=100, density=False)
    print('Histogram saved to ', 'figs/'+key+'.png')
    plt.savefig('figs/'+key+'.png')
    plt.yscale('log')
    plt.savefig('figs/'+key+'_log.png')
    plt.clf()

