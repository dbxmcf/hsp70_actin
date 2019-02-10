#!/usr/bin/env python
import numpy as np
import pandas as pd
import h5py

def rebuild_triangle(arr, st_loc, mtx_info):
    st = st_loc[0]
    #loc = st_loc[1]
    chunk_st_a = st_loc[2]
    #print(chunk_st)
    chunk_st_b = st_loc[3]
    chunk_ct_a = st_loc[4]
    #print(chunk_ct)
    chunk_ct_b = st_loc[5]
    total_lines = mtx_info[0]
    mpi_size = mtx_info[1]
    num_chunks = np.sqrt(2*mpi_size)
    if not num_chunks.is_integer:
        print("num_chunks is not integer")
    mat_wu = np.zeros((total_lines,total_lines))
    #print(mat_wu.shape)
    num_whole_blocks = int(num_chunks*(num_chunks-1)/2)

    arr_list = np.split(arr, st[1:])
    slc_nwb = slice(0,num_whole_blocks)
    #print(chunk_st_a[slc_nwb])
    for a,csta,cstb,ccta,cctb in zip(arr_list[slc_nwb],chunk_st_a[slc_nwb],chunk_st_b[slc_nwb],chunk_ct_a[slc_nwb],chunk_ct_b[slc_nwb]):
         mat_wu[csta:csta+ccta,cstb:cstb+cctb] = a.reshape(ccta,cctb)

    np.set_printoptions(edgeitems=30, linewidth=100000, formatter=dict(float=lambda x: "%.3f" % x))
    #print(mat_wu)

    #print("------")
    slc_ntri = slice(num_whole_blocks,mpi_size)
    #print(slc_ntri)
    for a,csta,cstb,ccta,cctb in zip(arr_list[slc_ntri],
                                     chunk_st_a[slc_ntri],chunk_st_b[slc_ntri],
                                     chunk_ct_a[slc_ntri],chunk_ct_b[slc_ntri]):
        #print(a,a.shape)
        #print(csta,cstb,ccta,cctb)
        sub_mat_a = mat_wu[csta:csta+ccta,csta:csta+ccta]

        segment_a = int(ccta*(ccta-1)/2)
        sub_mat_a[np.triu_indices(ccta,1)] = a[:segment_a]
        #print(sub_mat_a)
        #print(segment_a)
        #print(sub_mat_a[np.triu_indices(ccta,1)])
        #print(a[:segment_a+1])
        sub_mat_b = mat_wu[cstb:cstb+cctb,cstb:cstb+cctb]

        #segment_b = int(cctb*(cctb-1)/2)
        sub_mat_b[np.triu_indices(cctb,1)] = a[segment_a:]
        #print(sub_mat_b)

    i_lower = np.tril_indices(total_lines, -1)
    mat_wu[i_lower] = mat_wu.T[i_lower]
    print(mat_wu)
    #print(num_whole_blocks)
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
wu = np.array(f['wu'])[0]
#print("wu=",wu)
root_grp = f['/']
mtx_info = np.array(root_grp.attrs['MatrixInfo'])


rebuild_triangle(wu,start_loc,mtx_info)