#!/usr/bin/env
import numpy as np
from scipy import spatial
from numpy import linalg as LA
#a = np.array([1500,1501,1502,1503,1504,1505,1506,1507,1508,1509,1510,1511])
a = np.array([ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]) #,dtype=np.int32)
b = np.array([1600,1601,1602,1603,1604,1605,1606,1607,1608,1609,1610,1611]) #,dtype=np.int32)
#a = np.array([ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],dtype=np.int32)
#b = np.array([1600,1601,1602,1603,1604,1605,1606,1607,1608,1609,1610,1611],dtype=np.int32)
one_data_norm_a = 1.0/LA.norm(a)
one_data_norm_b = 1.0/LA.norm(b)
print(one_data_norm_a,one_data_norm_b)
print(a.dtype)
print(one_data_norm_a.dtype)
result = 1 - spatial.distance.cosine(a, b)
print("result=",result)
#print("adb=",a.dot(b))
#result1 = 1 - a.dot(b)*one_data_norm_a*one_data_norm_b
result1 = a.dot(b)*one_data_norm_a*one_data_norm_b
print("result1=",result1)
