#!/usr/bin/evn python

import numpy as np
import StringIO
import itertools
#import cartesian
from cartesian import *
from itertools import combinations

#f = "sample_hsp70_actin/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"
f = "test.csv"

#lines = np.loadtxt(f,usecols=(1,))
s = open(f).read().replace(';',',')
data = np.genfromtxt(StringIO.StringIO(s),delimiter=",")[:,1:]

#print(data)
print(data.shape)
print(data.dtype)

data_sum = np.sum(data,axis=1)
data_jac = data
data_jac[data>0]=1

#lst_a = np.arange(data.shape[0]).tolist()
#lst_b = np.arange(data.shape[0]).tolist()
#lst_a = np.arange(data.shape[0])
lst_a = np.arange(5)

print(list(combinations(lst_a,2)))
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


