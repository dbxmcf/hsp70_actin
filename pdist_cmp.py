#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import pandas as pd
import time
import numpy as np
#import StringIO
import itertools
from scipy import spatial
from numpy import linalg as LA
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
sample_name = "sample_a-b_mix_2"
#sample_name = "sample_protease_mix_1"
fname = sample_name + "/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"

start_time=time.time()

arrs = []
m_datatype = np.uint16

with open(fname) as fcsv:
    lines=fcsv.readlines()
    #for idx,line in enumerate(lines[n_st:]):
    for idx,line in enumerate(lines):
        l = list(line.split(';')[1].split(','))
        #l_arr = np.asarray(l[:-1]).astype(np.float) 
        l_arr = np.asarray(l[:-1],dtype=m_datatype)
        arrs.append(l_arr)
data = np.array(arrs,dtype=m_datatype)
print(data.shape)
print(data.dtype)

end_time=time.time()
total_time=((end_time)-(start_time))
print("Time taken for making matrix: {}".format(total_time))

#exit()

start_time=time.time()

cvt_type = np.float32
data_sum = np.sum(data,axis=1,dtype=cvt_type)
#data_jac = np.copy(data)
#data_jac[data_jac>0]=1
data_jac = data > 0

lst_a = np.arange(data.shape[0])

normal = np.zeros((data.shape[0],data.shape[0]))
generalised = np.zeros_like(normal)
sarika = np.zeros_like(normal)
wu = np.zeros_like(normal)
cosine = np.ones_like(normal)

lst_cmb = list(combinations(lst_a,2))
total_cmb = len(lst_cmb)
print("total_cmb=",total_cmb)
nitvl = min(total_cmb, 20)
itvl = total_cmb/nitvl
print("itvl=",itvl)

from scipy.spatial.distance import pdist
start_time=time.time()
jcd = pdist(data, metric='cosine')
cosine = 1.0-jcd
end_time=time.time()
total_time=((end_time)-(start_time))

csv_folder_name = sample_name+"_csv_"+m_datatype.__name__
pd.DataFrame(cosine).to_csv(csv_folder_name + "/pdist.csv")
print("Time taken for pdist jcd: {}".format(total_time))

#start_time=time.time()
#my_cosine=np.zeros(len(lst_cmb))
#for i, c in enumerate(lst_cmb):
#    idx_a, idx_b = c
#    a = data[idx_a]
#    b = data[idx_b]
#    result = 1 - spatial.distance.cosine(a.astype(cvt_type), b.astype(cvt_type))
#    my_cosine[i] = result
#pd.DataFrame(cosine).to_csv(csv_folder_name + "/my_pdist.csv")
#end_time=time.time()
#total_time=((end_time)-(start_time))
#print("Time taken for mypdist jcd: {}".format(total_time))

my_cosine=np.zeros(len(lst_cmb))
start_time=time.time()
one_data_norm = 1.0/LA.norm(data,axis=1)
for i, c in enumerate(lst_cmb):
    idx_a, idx_b = c
    a = data[idx_a]
    b = data[idx_b]
    one_an = one_data_norm[idx_a]
    one_bn = one_data_norm[idx_b]
    result = 1 - a.dot(b)*one_an*one_bn
    my_cosine[i] = result
pd.DataFrame(cosine).to_csv(csv_folder_name + "/my_opt_pdist.csv")
end_time=time.time()
total_time=((end_time)-(start_time))
print("Time taken for mypdist opt jcd: {}".format(total_time))



