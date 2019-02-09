#!/usr/bin/env python
import numpy as np
import pandas as pd
import h5py

def rebuild_triangle(arr, st_loc, mtx_info):
    st = st_loc[0]
    loc = st_loc[1]
    print(st)
    print(loc)
    total_lines = mtx_info[0]
    mpi_size = mtx_info[1]
    num_chunks = np.sqrt(2*mpi_size)
    if not num_chunks.is_integer:
        print("num_chunks is not integer")
    mat_wu = np.zeros((total_lines,total_lines))
    num_whole_blocks = int(num_chunks*(num_chunks-1)/2)
    print(num_whole_blocks)
    #for ():
    #    mat_wu[][] = 5



f = h5py.File('res_all.h5', 'r')
keys = list(f.keys())
print("keys=",keys)
start_loc = np.array(f['start_loc'])
#print("start_loc=",start_loc)
sarika = np.array(f['sarika'])
#print("sarika=",sarika)
normal = np.array(f['normal'])
#print("normal=",normal)
generalised = np.array(f['generalised'])
#print("generalised=",generalised)
cosine = np.array(f['cosine'])
#print("cosine=",cosine)
wu = np.array(f['wu'])
#print("wu=",wu)
root_grp = f['/']
mtx_info = np.array(root_grp.attrs['MatrixInfo'])


rebuild_triangle(wu,start_loc,mtx_info)