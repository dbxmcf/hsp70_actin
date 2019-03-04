#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import pandas as pd
import time
import numpy as np
import StringIO
import itertools
from scipy import spatial
#import cartesian
from cartesian import *
from itertools import combinations

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

start_time=time.time()

arrs = []
#m_datatype = np.uint8
m_datatype = np.uint16

with open(fname) as fcsv:
    lines=fcsv.readlines()
    #for idx,line in enumerate(lines[n_st:]):
    for idx,line in enumerate(lines):
        l = list(line.split(';')[1].split(','))
        #l_arr = np.asarray(l[:-1]).astype(np.float) 
        l_arr = np.asarray(l[:-1],dtype=m_datatype)
        arrs.append(l_arr)
        print('max=',l_arr.max(),'min=',l_arr.min(),'size=',l_arr.size)
        #print('max=',l_arr.max()) #,'min=',l_arr.min(),'size=',l_arr.size())
        #print('min=',l_arr.min())
        #print(
#data = np.array(arrs,dtype=m_datatype)
#print(data.shape)
#print(data.dtype)

end_time=time.time()
total_time=((end_time)-(start_time))
print("Time taken for making matrix: {}".format(total_time))

#exit()

start_time=time.time()
