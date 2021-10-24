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
dt = f1['oscw']
du = f1['weight']
dv = f1['ipnu']
ityp = np.array(ds[()])
osc = np.array(dt[()])
wgh = np.array(du[()])
neu = np.array(dv[()])
wo = osc*wgh

# for key in f1.keys():
# #    print(key)
#     ds = f1[key]
#     dat = np.array(ds[()])
# #    print(data)
#     for i in range(16):
#         for nu in [-16,-14,-12,12,14,16]:
#             data = dat[neu==nu]
#             w = wo[neu==nu]
#             typ = ityp[neu==nu]
#             # plt.hist(data[typ==i-1], weights=w[typ==i-1], bins=20, density=False, stacked=True, label=str(nu))
#             plt.hist(data[typ==i-1], bins=50, density=False, stacked=True, label=str(nu))
#         plt.legend()
#         plt.savefig('figs/'+key+'_'+str(i-1)+'.png')
#         print('Histogram saved to ', 'figs/'+key+'.png')
#         # plt.yscale('log')
#         # plt.savefig('figs/'+key+'_'+str(i-1)+'_log.png')
#         plt.clf()

ds = f1['dirnuZ']
true_cz = np.array(ds[()])
ds = f1['recodirZ']
reco_cz = np.array(ds[()])
for i in range(16):
    for nu in [-16,-14,-12,12,14,16]:
        cond = (ityp==i)*(np.absolute(true_cz)>0.8)*(neu==nu)
        plt.hist(true_cz[cond]-reco_cz[cond], bins=20, density=False)
        plt.title(str(i)+str(nu))
        # plt.hist(true_cz[ityp==i-1], bins=10, weights=wo[ityp==i-1], density=False)
        # plt.show()
        # plt.clf()
        # plt.hist(reco_cz[ityp==i-1], bins=10, weights=wo[ityp==i-1], density=False)
        plt.show()
        plt.clf()
