#!/usr/bin/env python

from __future__ import print_function
import pandas as pd
import time
import numpy as np
import StringIO
import itertools
from scipy import spatial
#import cartesian
from cartesian import *
from itertools import combinations

f = "sample_hsp70_actin/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"
#f = "test.csv"

#lines = np.loadtxt(f,usecols=(1,))
s = open(f).read().replace(';',',')
data = np.genfromtxt(StringIO.StringIO(s),delimiter=",")[:,1:-1]

start_time=time.time()
#print(data)
print(data.shape)
print(data.dtype)

data_sum = np.sum(data,axis=1)
data_jac = np.copy(data)
data_jac[data_jac>0]=1

#lst_a = np.arange(data.shape[0]).tolist()
#lst_b = np.arange(data.shape[0]).tolist()
lst_a = np.arange(data.shape[0])
#lst_a = np.arange(5)

lst_cmb = list(combinations(lst_a,2))
print(len(lst_cmb))

normal = np.zeros((data.shape[0],data.shape[0]))
generalised = np.zeros_like(normal)
sarika = np.zeros_like(normal)
wu = np.zeros_like(normal)
cosine = np.zeros_like(normal)

for c in lst_cmb:
    idx_a, idx_b = c
    #print("a=",a,"b=",b)
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
        
pd.DataFrame(normal).to_csv("csv/normal.csv")
pd.DataFrame(generalised).to_csv("csv/generalised.csv")
pd.DataFrame(sarika).to_csv("csv/sarika1.csv")
pd.DataFrame(wu).to_csv("csv/wu.csv")
pd.DataFrame(cosine).to_csv("csv/cosine.csv")

end_time=time.time()
total_time=((end_time)-(start_time))
print("Time taken for writing to files: {}".format(total_time))


    #    normal.append(str(dist_jac)) 
    #    generalised.append(str(dist_gen_jac))                      
    #    sarika.append(str(dist_sarika))
    #    wu.append(str(dist_wu))
    #    cosine.append(str(result*100))
    #return (name, normal, generalised, sarika, wu, cosine)
#lst_b = np.arange(data.shape[0])
#
#print(lst_a)
#print(lst_b)
#
#permutation_list = cartesian((lst_a, lst_b))
#print(permutation_list)
#non_zeros = 

#            i=list(i.split(';')[1].split(','))
#            a = np.asarray(i[:-1]).astype(np.float)  
#            a_sum = np.sum(a)
#            a_jac = np.copy(a)
#            a_jac[a_jac>0] = 1
#            
#            for j in self.lines:
#                j = list(j.split(';')[1].split(','))
#                b = np.asarray(j[:-1]).astype(np.float)
#                non_zeros = (a >0) & (b > 0)
#                summed_array = a + b
#                b_sum = np.sum(b)
#                b_jac = np.copy(b)
#                b_jac[b_jac>0] = 1
#                
#                numerator_jac = np.sum(np.minimum(a_jac,b_jac))
#                denomenator_jac = np.sum(np.maximum(a_jac,b_jac))
#                numerator_gen_jac =np.sum(np.minimum(a,b))
#                denomenator_gen_jac =np.sum(np.maximum(a,b))
#                num_sim = np.sum(summed_array[non_zeros])
#                result = 1 - spatial.distance.cosine(a, b)


