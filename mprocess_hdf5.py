import os
import numpy as np
from multiprocessing import Pool
import h5py
def distance_block(idx):
    """ Read a 100-element coordinates block, compute distances, and write
    back out again to a process-specific file.
    """
    with h5py.File('coords.hdf5','r') as f:
        data = f['coords'][idx:idx+100]
    result = np.sqrt(np.sum(data**2, axis=1))
    with h5py.File('result_index_%d.hdf5'%idx, 'w') as f:
        f['result'] = result

# Create out pool and carry out the computation
p = Pool(4)
p.map(distance_block, np.arange(0, 1000, 100))
with h5py.File('coords.hdf5') as f:
    dset = f.create_dataset('distances', (1000,), dtype='f4')
    # Loop over our 100-element "chunks" and merge the data into coords.hdf5
    for idx in xrange(0, 1000, 100):
        filename = 'result_index_%d.hdf5'%idx
        with h5py.File(filename, 'r') as f2:
            data = f2['result'][...]
        dset[idx:idx+100] = data
        os.unlink(filename) # no longer needed