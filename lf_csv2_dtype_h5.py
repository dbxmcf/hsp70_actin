#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import pandas as pd
import time
import numpy as np
from numpy import linalg as LA
#import StringIO
import itertools
from scipy import spatial
#import cartesian
#from cartesian import *
from itertools import combinations
import h5py

#orig_stdout = sys.stdout
#f = open('log.txt', 'w')
#sys.stdout = f

#fname = "test.csv"
#fname = "sample_hsp70_actin/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"
#fname = sample_name + "/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"
#fname = "sample_protease_mix_1/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"

#sample_name = "sample_hsp70_actin"
#sample_name = "sample_a-b_mix_2"
sample_name = "sample_protease_mix_1"
fname = sample_name + "/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"

#sample_name = "hdf5t"
#fname = sample_name + "/ta.csv"

start_time=time.time()

arrs = []
m_datatype = np.uint16
#m_datatype = np.int64

with open(fname) as fcsv:
    lines=fcsv.readlines()
    #for idx,line in enumerate(lines[n_st:]):
    for idx,line in enumerate(lines):
        l = list(line.split(';')[1].split(','))
        #l_arr = np.asarray(l[:-1]).astype(np.float) 
        l_arr = np.asarray(l[:-1],dtype=m_datatype)
        #l_arr = np.asarray(l[:],dtype=m_datatype)
        arrs.append(l_arr)
data = np.array(arrs,dtype=m_datatype)
#print(data)
print(data.shape)
print(data.dtype)
print("min val=",data.min())
print("max val=",data.max())

end_time=time.time()
total_time=((end_time)-(start_time))
print("Time taken for making matrix: {}".format(total_time))

h5_filename = "hdf5t/" + sample_name + "_dtype.h5"
#h5_filename = "hdf5t/" + sample_name + ".h5"
h5f = h5py.File(h5_filename, 'w')
#h5f.create_dataset('Data1', data=data, dtype='int16')
#h5f.create_dataset('Data1', data=data, dtype='int16')
h5f.create_dataset('Data1', data=data, dtype=m_datatype)
h5f.close()
print(h5_filename + " file created")
#exit()

