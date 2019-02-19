#!/usr/bin/env python
import sys
import numpy as np
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
    mat_ret = np.zeros((total_lines,total_lines))
    #print(mat_ret.shape)
    num_whole_blocks = int(num_chunks*(num_chunks-1)/2)

    arr_list = np.split(arr, st[1:])
    slc_nwb = slice(0,num_whole_blocks)
    #print(chunk_st_a[slc_nwb])
    for a,csta,cstb,ccta,cctb in zip(arr_list[slc_nwb],chunk_st_a[slc_nwb],chunk_st_b[slc_nwb],chunk_ct_a[slc_nwb],chunk_ct_b[slc_nwb]):
         mat_ret[csta:csta+ccta,cstb:cstb+cctb] = a.reshape(ccta,cctb)

    np.set_printoptions(edgeitems=30, linewidth=100000, formatter=dict(float=lambda x: "%.3f" % x))
    #print(mat_ret)

    #print("------")
    slc_ntri = slice(num_whole_blocks,mpi_size)
    #print(slc_ntri)
    for a,csta,cstb,ccta,cctb in zip(arr_list[slc_ntri],
                                     chunk_st_a[slc_ntri],chunk_st_b[slc_ntri],
                                     chunk_ct_a[slc_ntri],chunk_ct_b[slc_ntri]):
        #print(a,a.shape)
        #print(csta,cstb,ccta,cctb)
        sub_mat_a = mat_ret[csta:csta+ccta,csta:csta+ccta]

        segment_a = int(ccta*(ccta-1)/2)
        sub_mat_a[np.triu_indices(ccta,1)] = a[:segment_a]
        #print(sub_mat_a)
        #print(segment_a)
        #print(sub_mat_a[np.triu_indices(ccta,1)])
        #print(a[:segment_a+1])
        sub_mat_b = mat_ret[cstb:cstb+cctb,cstb:cstb+cctb]

        #segment_b = int(cctb*(cctb-1)/2)
        sub_mat_b[np.triu_indices(cctb,1)] = a[segment_a:]
        #print(sub_mat_b)

    i_lower = np.tril_indices(total_lines, -1)
    mat_ret[i_lower] = mat_ret.T[i_lower]
    return mat_ret
    #print(mat_wu)
    #print(num_whole_blocks)
    #for ():
    #    mat_wu[][] = 5



f = h5py.File('res_all.h5', 'r')
keys = list(f.keys())
print("keys=",keys)
start_loc = np.array(f['start_loc'])
#print("start_loc=",start_loc)
sarika = np.array(f['sarika'])[0]
#print("sarika=",sarika)
normal = np.array(f['normal'])[0]
#print("normal=",normal)
generalised = np.array(f['generalised'])[0]
#print("generalised=",generalised)
cosine = np.array(f['cosine'])[0]
#print("cosine=",cosine)
wu = np.array(f['wu'])[0]
#print("wu=",wu)
root_grp = f['/']
mtx_info = np.array(root_grp.attrs['MatrixInfo'])


mat_normal_h5 = rebuild_triangle(normal,start_loc,mtx_info)
mat_wu_h5 = rebuild_triangle(wu,start_loc,mtx_info)
mat_generalised_h5 = rebuild_triangle(generalised,start_loc,mtx_info)
mat_sarika_h5 = rebuild_triangle(sarika,start_loc,mtx_info)
mat_cosine_h5 = rebuild_triangle(cosine,start_loc,mtx_info)

sample_name = sys.argv[1]
py_path_name = "../" + sample_name +"_csv_uint16/"
#py_path_name = "../hdf5t_csv_uint16/"
#py_path_name = "../sample_hsp70_actin_csv_uint16/"
#py_path_name = "../sample_a-b_mix_2_csv_uint16/"
#py_path_name = "../sample_protease_mix_1_csv_uint16/"
# for key in keys:
mat_normal_py = np.loadtxt(py_path_name + "normal.csv",delimiter=",")
mat_tol = 0.002
ret = np.allclose(mat_normal_h5, mat_normal_py, atol=mat_tol)
if (ret):
    print("normal is ok!")
else:
    print("normal is not ok!")

# for key in keys:
mat_sarika_py = np.loadtxt(py_path_name + "sarika.csv",delimiter=",")
ret = np.allclose(mat_sarika_h5, mat_sarika_py, atol=mat_tol)
if (ret):
    print("sarika is ok!")
else:
    print("sarika is not ok!")
    #np.savetxt('sarika.csv',fmt="%7.3f")
    np.savetxt("sarika.csv", mat_sarika_h5, delimiter=",",fmt="%7.3f")

mat_generalised_py = np.loadtxt(py_path_name + "generalised.csv",delimiter=",")
ret = np.allclose(mat_generalised_h5, mat_generalised_py, atol=mat_tol)
if (ret):
    print("generalised is ok!")
else:
    print("generalised is not ok!")

mat_wu_py = np.loadtxt(py_path_name + "wu.csv",delimiter=",")
ret = np.allclose(mat_wu_h5, mat_wu_py, atol=mat_tol)
if (ret):
    print("wu is ok!")
else:
    print("wu is not ok!")

#mat_cosine_py = np.loadtxt(py_path_name + "cosine.csv",delimiter=",")
#ret = np.allclose(mat_cosine_h5, mat_cosine_py, atol=mat_tol)
#if (ret):
#    print("cosine is ok!")
#else:
#    print("cosine is not ok!")
#print(mat_normal_py)
#np.set_printoptions(edgeitems=30, linewidth=100000, formatter=dict(float=lambda x: "%7.3f" % x))
#print(mat_cosine)
