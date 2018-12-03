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
sample_name = "sample_a-b_mix_2"
#sample_name = "sample_protease_mix_1"
fname = sample_name + "/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"

start_time=time.time()

arrs = []
m_datatype = np.float32

with open(fname) as fcsv:
    lines=fcsv.readlines()
    n_lines = len(lines)
    n_div = 2
    n_st = n_lines/n_div
    #for idx,line in enumerate(lines[n_st:]):
    for idx,line in enumerate(lines):
        l = list(line.split(';')[1].split(','))
        #l_arr = np.asarray(l[:-1]).astype(np.float) 
        l_arr = np.asarray(l[:-1],dtype=m_datatype)
        arrs.append(l_arr)
data = np.array(arrs,dtype=m_datatype)
del arrs[:]
print(data.shape)
print(data.dtype)

end_time=time.time()
total_time=((end_time)-(start_time))
print("Time taken for making matrix: {}".format(total_time))

#exit()

start_time=time.time()

data_sum = np.sum(data,axis=1)
data_jac = np.copy(data)
data_jac[data_jac>0]=1

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

for i, c in enumerate(lst_cmb):
    idx_a, idx_b = c
    a = data[idx_a]
    a_sum = data_sum[idx_a]
    a_jac = data_jac[idx_a]
    b = data[idx_b]
    b_sum = data_sum[idx_b]
    b_jac = data_jac[idx_b]

    non_zeros = (a >0) & (b > 0)
    summed_array = a + b

    numerator_jac = np.sum(np.minimum(a_jac,b_jac))
    denomenator_jac = np.sum(np.maximum(a_jac,b_jac))
    numerator_gen_jac =np.sum(np.minimum(a,b))
    denomenator_gen_jac =np.sum(np.maximum(a,b))
    num_sim = np.sum(summed_array[non_zeros])
    result = 1 - spatial.distance.cosine(a, b)

    if (denomenator_jac == 0):
        print('There is something wrong. Denominator is Zero! ', idx_a, idx_b, numerator_jac, denomenator_jac)
    else:
        dist_gen_jac=1.0-(float(numerator_gen_jac)/float(denomenator_gen_jac))                    
        dist_jac=1.0-(float(numerator_jac)/float(denomenator_jac))

        denomenator_wu = min(denomenator_gen_jac,max(a_sum,b_sum) )
        dist_wu = 1.0-(float(numerator_gen_jac)/float(denomenator_wu))
        
        numerator_sarika = num_sim
        denomenator_sarika = a_sum+b_sum
        dist_sarika = 1.0-(float(numerator_sarika)/float(denomenator_sarika))

    normal[idx_a,idx_b] = dist_jac
    normal[idx_b,idx_a] = dist_jac
    generalised[idx_a,idx_b] = dist_gen_jac
    generalised[idx_b,idx_a] = dist_gen_jac
    sarika[idx_a,idx_b] = dist_sarika
    sarika[idx_b,idx_a] = dist_sarika
    wu[idx_a,idx_b] = dist_wu
    wu[idx_b,idx_a] = dist_wu
    cosine[idx_a,idx_b] = result*100
    cosine[idx_b,idx_a] = result*100
    if (i % itvl) == 0:
        print("iter:\t",i,"\ttime at {}".format(time.time()-start_time))
        
csv_folder_name = sample_name+"_csv_"+m_datatype.__name__
if not os.path.exists(csv_folder_name):
    os.mkdir(csv_folder_name)
pd.DataFrame(normal).to_csv(csv_folder_name+"/normal.csv")
pd.DataFrame(generalised).to_csv(csv_folder_name+"/generalised.csv")
pd.DataFrame(sarika).to_csv(csv_folder_name+"/sarika1.csv")
pd.DataFrame(wu).to_csv(csv_folder_name + "/wu.csv")
pd.DataFrame(cosine).to_csv(csv_folder_name + "/cosine.csv")

end_time=time.time()
total_time=((end_time)-(start_time))
print("Time taken for Jaccard: {}".format(total_time))

#sys.stdout = orig_stdout
#f.close()

