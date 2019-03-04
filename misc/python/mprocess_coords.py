import numpy as np
from multiprocessing import Pool
import h5py
def distance(arr):
    """ Compute distance from origin to the point (arr is a shape-(2,) array)
    """
    return np.sqrt(np.sum(arr**2))
# Load data and close the input file
with h5py.File('coords.hdf5', 'r') as f:
    data = f['coords'][...]
# Create a 4-process pool
p = Pool(4)
# Carry out parallel computation
result = np.array(p.map(distance, data))
# Write the result into a new dataset in the file
with h5py.File('coords.hdf5') as f:
    f['distances'] = result