import h5py
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('fname', type=str, nargs='?', default='NULL')
parser.add_argument('variable', type=str, nargs='?', default='ALL')
args = parser.parse_args()
infile = args.fname
var = args.variable
print(var)

if infile == 'NULL':
    sys.exit('Please introduce your HDF5 file')

if var != 'ALL':
    f1 = h5py.File(infile, 'r')
    ds = f1[var]
    data = np.array(ds[()])
    typ = np.array(f1['itype'][()])
    for i in range(16):
        plt.hist(data[typ==i], bins=50, density=False, stacked=True)
        plt.show()
        plt.clf()
else:
    for key in f1.keys():
    #    print(key)
        ds = f1[key]
        dat = np.array(ds[()])
    #    print(data)
        for i in range(16):
            for nu in [-16,-14,-12,12,14,16]:
                data = dat[neu==nu]
                w = wo[neu==nu]
                typ = ityp[neu==nu]
                # plt.hist(data[typ==i-1], weights=w[typ==i-1], bins=20, density=False, stacked=True, label=str(nu))
                plt.hist(data[typ==i-1], bins=50, density=False, stacked=True, label=str(nu))
            plt.legend()
            plt.savefig('../figs/'+key+'_'+str(i-1)+'.png')
            print('Histogram saved to ', 'figs/'+key+'.png')
            # plt.yscale('log')
            # plt.savefig('figs/'+key+'_'+str(i-1)+'_log.png')
            plt.clf()

