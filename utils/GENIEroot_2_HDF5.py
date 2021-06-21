import uproot as upt
import h5py
import numpy as np
import argparse
import sys
#import awkward as ak

parser = argparse.ArgumentParser()
parser.add_argument('fname', type=str, nargs='?', default='NULL')
args = parser.parse_args()
infile = args.fname

if infile == 'NULL':
    sys.exit('Hey! You forgot to introduce your GENIE file')

f = upt.open(infile)
keep_columns = f['gst'].keys()
print(keep_columns)

dt = h5py.special_dtype(vlen=np.float32)

with h5py.File(infile+'.hdf5', 'w') as hf:
    for i,br in enumerate(keep_columns):
        dummy = f.get('gst')[br].array(interpretation=None, entry_start=None, entry_stop=None, decompression_executor=None, interpretation_executor=None, array_cache='inherit', library='np')
#        pd_dummy = ak.to_pandas(dummy)
#        pd_dummy.to_hdf(h5File, br);
#        print(dummy)
        dummy = np.array(dummy)
        if dummy.dtype==np.dtype('object'):
            # Special treatment for variable length datasets
            dt = h5py.special_dtype(vlen=np.float64)
            hf.create_dataset(br, data=dummy, dtype=dt, compression='gzip')
#            hf.create_dataset(br, data=dummy, dtype=dt)
        else:
            hf.create_dataset(br, data=dummy, compression='gzip')
#            hf.create_dataset(br, data=dummy)


