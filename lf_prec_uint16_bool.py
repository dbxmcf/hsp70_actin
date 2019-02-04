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
#sample_name = "sample_protease_mix_1"
#fname = sample_name + "/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"

sample_name = "hdf5t"
fname = sample_name + "/ta.csv"

start_time=time.time()

arrs = []
m_datatype = np.uint16

with open(fname) as fcsv:
    lines=fcsv.readlines()
    #for idx,line in enumerate(lines[n_st:]):
    for idx,line in enumerate(lines):
        l = list(line.split(';')[1].split(','))
        #l_arr = np.asarray(l[:-1]).astype(np.float) 
        l_arr = np.asarray(l[:],dtype=m_datatype)
        arrs.append(l_arr)
data = np.array(arrs,dtype=m_datatype)
print(data)
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

data1=data.astype(np.float64)
one_data_norm = 1.0/LA.norm(data1,axis=1)

print("one_data_norm:",one_data_norm)
for i, c in enumerate(lst_cmb):
    idx_a, idx_b = c
    a = data[idx_a]
    a_sum = data_sum[idx_a]
    a_jac = data_jac[idx_a]
    b = data[idx_b]
    b_sum = data_sum[idx_b]
    b_jac = data_jac[idx_b]

    #non_zeros = (a >0) & (b > 0)
    non_zeros = a_jac & b_jac
    summed_array = a + b

    #numerator_jac = np.sum(np.minimum(a_jac,b_jac))
    #denomenator_jac = np.sum(np.maximum(a_jac,b_jac))
    #numerator_jac = np.sum(non_zeros)
    #denomenator_jac = np.sum(a_jac | b_jac)
    numerator_jac = np.count_nonzero(non_zeros)
    denomenator_jac = np.count_nonzero(a_jac | b_jac)
    numerator_gen_jac =np.sum(np.minimum(a,b))
    denomenator_gen_jac =np.sum(np.maximum(a,b))
    num_sim = np.sum(summed_array[non_zeros])
    #result = 1 - spatial.distance.cosine(a.astype(cvt_type), b.astype(cvt_type))
    #result = 1 - spatial.distance.cosine(a, b)
    one_an = one_data_norm[idx_a]
    one_bn = one_data_norm[idx_b]
    #print(a)
    #print(b)
    a1 = a.astype(np.int)
    b1 = b.astype(np.int)
    adotb = a1.dot(b1)
    #print(idx_a, idx_b, adotb)
    one_an1 = one_an.astype(np.float64)
    one_bn1 = one_bn.astype(np.float64)
    a1dotb1 = a1.dot(b1)
    res1 = 1 - a1.dot(b1)*one_an1*one_bn1
    print(idx_a, idx_b, one_an1, one_bn1, a1dotb1, res1)
    #result = 1 - a1.dot(b1)*one_an1*one_bn1
    result = a1.dot(b1)*one_an1*one_bn1
    #result = 1 - spatial.distance.cosine(a.astype(cvt_type), b.astype(cvt_type))

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
    #print("normal=",dist_jac)
    #print("generalised=",dist_gen_jac)
    #print("sarika=",dist_sarika)
    #print("wu=",dist_wu)
    #print("cosine=",result*100)
    #if (i % itvl) == 0:
    #    print("iter:\t",i,"\ttime at {}".format(time.time()-start_time))
        
csv_folder_name = sample_name+"_csv_"+m_datatype.__name__
if not os.path.exists(csv_folder_name):
    os.mkdir(csv_folder_name)

csv_fmt = '%7.3f'
np.savetxt(csv_folder_name+"/normal.csv", normal, delimiter=",",fmt=csv_fmt)
np.savetxt(csv_folder_name+"/generalised.csv", generalised, delimiter=",",fmt=csv_fmt)
np.savetxt(csv_folder_name+"/sarika.csv", sarika, delimiter=",",fmt=csv_fmt)
np.savetxt(csv_folder_name+"/wu.csv", wu, delimiter=",",fmt=csv_fmt)
np.savetxt(csv_folder_name+"/cosine.csv", cosine, delimiter=",",fmt=csv_fmt)
#pd.DataFrame(normal).to_csv(csv_folder_name+"/normal.csv")
#pd.DataFrame(generalised).to_csv(csv_folder_name+"/generalised.csv")
#pd.DataFrame(sarika).to_csv(csv_folder_name+"/sarika1.csv")
#pd.DataFrame(wu).to_csv(csv_folder_name + "/wu.csv")
#pd.DataFrame(cosine).to_csv(csv_folder_name + "/cosine.csv")

end_time=time.time()
total_time=((end_time)-(start_time))
print("Time taken for Jaccard: {}".format(total_time))

#sys.stdout = orig_stdout
#f.close()

