#!/usr/bin/env python

from __future__ import print_function
from mpi4py import MPI
import sys
import os
import pandas as pd
import time
import numpy as np
#import itertools
from scipy import spatial
from itertools import combinations
from operator import itemgetter

def chunk_list(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

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
sample_name = "tmpi"
fname = sample_name + "/test.csv"

start_time=time.time()

arrs = []
m_datatype = np.float32

comm = MPI.COMM_WORLD
nprocs = comm.Get_size()
rank = comm.Get_rank()

# divide the entire list into n parts
n_parts = int(np.sqrt(2.0*nprocs))
#print("n_parts=", n_parts)

# off-diagnal combinations
offdiag_cmbs = list(combinations(range(n_parts),2))
#print(offdiag_cmbs)
# diagnal combinations
diag_cmbs = [(i, i) for i in range(n_parts)]
#print(diag_cmbs)
diag_cmbs_fold = []
n_half_parts = int(n_parts/2)
for i in range(n_half_parts):
    diag_cmbs_fold.append((diag_cmbs[i],diag_cmbs[-1-i]))
#print(diag_cmbs_fold)

parts_cmbs = offdiag_cmbs + diag_cmbs_fold
#print(parts_cmbs)
parts_cmb_rank = parts_cmbs[rank]
#if rank > int(n_parts*(n_parts-1)/2-1):
#print("rank=",rank,", parts_cmb_rank=", parts_cmb_rank)

with open(fname) as fcsv:
    lines=fcsv.readlines()
    n_lines = len(lines)

line_numbers = list(range(n_lines))
chunks = chunk_list(line_numbers, n_parts)

rank_div = int(n_parts*(n_parts-1)/2-1)

if rank <= rank_div:
    #print("rank=",rank,", parts_cmb_rank=", parts_cmb_rank)
    parts_lines = chunks[parts_cmb_rank[0]] + chunks[parts_cmb_rank[1]]
    parts_line_cmbs = list(combinations(parts_lines,2))
    #print("rank=",rank,", parts_cmb_rank=", parts_cmb_rank, ", parts_lines=", parts_lines,", parts_line_cmbs=", parts_line_cmbs) 
else:
    parts_lines = chunks[parts_cmb_rank[0][0]] + chunks[parts_cmb_rank[1][0]]
    parts_line_cmbs = (list(combinations(chunks[parts_cmb_rank[0][0]],2))
                       + list(combinations(chunks[parts_cmb_rank[1][0]],2)))
    #print("rank=",rank,", parts_cmb_rank=", parts_cmb_rank, ", parts_lines=", parts_lines,", parts_line_cmbs=", parts_line_cmbs) 

rank_lines = itemgetter(*parts_lines)(lines)

data = {}
data_sum = {}
data_jac = {}
for pl, line in zip(parts_lines,rank_lines):
    l = list(line.split(';')[1].split(','))
    #l_arr = np.asarray(l[:-1]).astype(np.float) 
    data[pl] = np.asarray(l[:-1],dtype=m_datatype)
    #arrs.append(l_arr)
    data_sum[pl] = np.sum(data[pl])
    data_jac[pl] = np.copy(data[pl])
    data_jac[pl][data_jac[pl]>0]=1

#data = np.array(arrs,dtype=m_datatype)

#print(data.shape)
#print(data.dtype)

end_time=time.time()
total_time=((end_time)-(start_time))
#print("Time taken for making matrix: {}".format(total_time))

#exit()

if rank == 0:
    start_time=time.time()


#for l, d in zip(parts_lines, data):
    #print("rank=",rank,l,d)
#    data_sum[l] = np.sum(d)
#    data_jac[l] = np.copy(d)
#    data_jac[l][data_jac[l]>0]=1

#print(rank, data_sum, data_jac)
#lst_a = np.arange(data.shape[0])
#exit()

#normal = np.zeros((len(parts_lines),len(parts_lines)))
#generalised = np.zeros_like(normal)
#sarika = np.zeros_like(normal)
#wu = np.zeros_like(normal)
#cosine = np.ones_like(normal)


total_cmbs_in_rank = len(parts_line_cmbs)
print("total_cmb_in_rank[", rank, "]=",total_cmbs_in_rank)
#nitvl = min(total_cmb, 20)
#itvl = total_cmb/nitvl
#print("itvl=",itvl)
rank_idx = np.zeros((total_cmbs_in_rank,2), dtype = int)
normal = np.zeros(total_cmbs_in_rank, dtype=m_datatype)
generalised = np.zeros_like(normal, dtype=m_datatype)
sarika = np.zeros_like(normal, dtype=m_datatype)
wu = np.zeros_like(normal, dtype=m_datatype)
cosine = np.zeros_like(normal, dtype=m_datatype)

#if rank > int(n_parts*(n_parts-1)/2-1):
#    exit()
for idx, c in enumerate(parts_line_cmbs):
    idx_a, idx_b = c
    #print(idx_a, idx_b)
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
    
    #normal.append((idx_a, idx_b, dist_jac))
    #generalised.append((idx_a, idx_b, dist_gen_jac))
    #sarika.append((idx_a, idx_b, dist_sarika))
    #wu.append((idx_a, idx_b, dist_wu))
    #cosine.append((idx_a, idx_b, result*100))
    
    #rank_idx[idx,0], rank_idx[idx,1] = idx_a, idx_b
    normal[idx] = dist_jac
    generalised[idx] = dist_gen_jac
    sarika[idx] = dist_sarika
    wu[idx] = dist_wu
    cosine[idx] = result*100
    
    #normal[idx_a,idx_b] = dist_jac
    #normal[idx_b,idx_a] = dist_jac
    #generalised[idx_a,idx_b] = dist_gen_jac
    #generalised[idx_b,idx_a] = dist_gen_jac
    #sarika[idx_a,idx_b] = dist_sarika
    #sarika[idx_b,idx_a] = dist_sarika
    #wu[idx_a,idx_b] = dist_wu
    #wu[idx_b,idx_a] = dist_wu
    #cosine[idx_a,idx_b] = result*100
    #cosine[idx_b,idx_a] = result*100
    #if (i % itvl) == 0:
    #    print("itvl:\t",i,"\ttime at {}".format(time.time()-start_time))

# need to do gather for data here

#exit()
        
#if rank == 0:
csv_folder_name = sample_name+"_csv_"+m_datatype.__name__
if not os.path.exists(csv_folder_name):
    os.mkdir(csv_folder_name)
pd.DataFrame(normal).to_csv(csv_folder_name+"/normal_" + str(rank) + ".csv")
pd.DataFrame(generalised).to_csv(csv_folder_name+"/generalised_" + str(rank) + ".csv")
pd.DataFrame(sarika).to_csv(csv_folder_name+"/sarika_" + str(rank) + ".csv")
pd.DataFrame(wu).to_csv(csv_folder_name + "/wu_" + str(rank) + ".csv")
pd.DataFrame(cosine).to_csv(csv_folder_name + "/cosine_" + str(rank) + ".csv")

if rank == 0:
    end_time=time.time()
    total_time=((end_time)-(start_time))
    print("Time taken for Jaccard: {}".format(total_time))

#sys.stdout = orig_stdout
#f.close()

