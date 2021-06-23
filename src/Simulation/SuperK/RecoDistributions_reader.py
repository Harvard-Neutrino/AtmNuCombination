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
    sys.exit('Please introduce your SK distributions HDF5 file')

f1 = h5py.File(infile, 'r')

for key in f1.keys():
#    print(key)
    ds = f1[key]['values']
    d = np.array(ds[()])
    ds = f1[key]['weights']
    w = np.array(ds[()])
#    print(data)
#    if data.dtype == np.dtype('object'):
#        plt.hist(np.hstack(data.flatten()))
#    elif data.dtype != np.dtype('bool'):
    plt.hist(d, bins=w.size, weights=w, density=True)
    print('Histogram saved to ', 'figs/'+key+'.png')
    plt.savefig('figs/'+key+'.png')
    plt.clf()

