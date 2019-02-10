#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import pandas as pd

#sample_name = "sample_hsp70_actin"
#sample_name = "sample_a-b_mix_2"
sample_name = "sample_protease_mix_1"
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

#complib = 'zlib'
complib = 'bzip2'
filename = sample_name + '_' + complib + '.h5'

#df = pd.DataFrame(np.arange(10).reshape((5,2)), columns=['A', 'B'])
#df = pd.DataFrame(np.arange(10).reshape((5,2)))
df = pd.DataFrame(data)
#print(df)
#    A  B
# 0  0  1
# 1  2  3
# 2  4  5
# 3  6  7
# 4  8  9

# Save to HDF5
#df.to_hdf(filename, 'data', mode='w', format='table')
df.to_hdf(filename, 'data', mode='w', complib=complib, complevel=9) #, format='table')
print('done with: ', filename)
del df    # allow df to be garbage collected
# Append more data
#df2 = pd.DataFrame(np.arange(10).reshape((5,2))*10, columns=['A', 'B'])
#df2.to_hdf(filename, 'data', append=True)

#print(pd.read_hdf(filename, 'data'))
