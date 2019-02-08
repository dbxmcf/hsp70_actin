#!/usr/bin/env python
import numpy as np
import pandas as pd
import h5py

def rebuild_triangle(arr, st_loc):
    st = st_loc[0]
    loc = st_loc[1]
    print(st)
    print(loc)

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

rebuild_triangle(wu,start_loc)