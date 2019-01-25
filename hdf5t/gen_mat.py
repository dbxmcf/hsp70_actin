#!/usr/bin/env python

#!/usr/bin/env python
import numpy as np
import pandas as pd
import h5py

val_range = 16
n = 17
m = 24
test_arr = np.random.randint(val_range, size=(n, m),dtype=np.uint16)
h5f = h5py.File('ta.h5', 'w')
h5f.create_dataset('test_arr_name', data=test_arr)
h5f.close()

print(test_arr)

exit()

sample_name = "sample_hsp70_actin"
#sample_name = "sample_a-b_mix_2"
#sample_name = "sample_protease_mix_1"
fname = sample_name + "/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"

m_datatype = np.uint16

arrs=[]
with open(fname) as fcsv:
    lines=fcsv.readlines()
    #for idx,line in enumerate(lines[n_st:]):
    for idx,line in enumerate(lines):
        l = list(line.split(';')[1].split(','))
        #l_arr = np.asarray(l[:-1]).astype(np.float) 
        l_arr = np.asarray(l[:-1],dtype=m_datatype)
        arrs.append(l_arr)
data = np.array(arrs,dtype=m_datatype)

filename = sample_name + '.h5'

np.save(sample_name + '.npy',data)
np.savetxt(sample_name + '.txt',data)

h5f = h5py.File(sample_name + '.h5', 'w')
h5f.create_dataset('ds1', data=data)
h5f.close()

#df = pd.DataFrame(np.arange(10).reshape((5,2)), columns=['A', 'B'])
#df = pd.DataFrame(np.arange(10).reshape((5,2)))
#df = pd.DataFrame(data)
#print(df)
#    A  B
# 0  0  1h5ls

# 1  2  3
# 2  4  5
# 3  6  7
# 4  8  9

# Save to HDF5
#df.to_hdf(filename, 'data', mode='w', format='table')
#df.to_hdf(filename, 'data', mode='w', complib=clib, complevel=9) #, format='table')
#del df    # allow df to be garbage collected
print('done with: ', filename)
# Append more data
#df2 = pd.DataFrame(np.arange(10).reshape((5,2))*10, columns=['A', 'B'])
#df2.to_hdf(filename, 'data', append=True)

#print(pd.read_hdf(filename, 'data'))
