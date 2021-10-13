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
ds = f1['itype']
typ = np.array(ds[()])

for key in f1.keys():
#    print(key)
    ds = f1[key]
    data = np.array(ds[()])
#    print(data)
    for i in range(20):
        plt.hist(data[(typ==i-1)], bins=100, density=False)
        print('Histogram saved to ', 'figs/'+key+'.png')
        plt.savefig('figs/'+key+'_'+str(i-1)+'.png')
        plt.yscale('log')
        plt.savefig('figs/'+key+'_'+str(i-1)+'_log.png')
        plt.clf()

