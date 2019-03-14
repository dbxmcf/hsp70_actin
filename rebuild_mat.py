#!/usr/bin/env python
import sys
import numpy as np
import h5py
import errno    
import os
import argparse
import json

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def save_mat(arr, fmt, lst_ptn_names, out_csv_fname):
    n_lines = arr.shape[0]
    with open(out_csv_fname, 'w') as f:
        for i in np.arange(0,n_lines):
            #str_arr_i = ','.join(map(str, arr[i]))
            #str_arr_i = ','.join(map("{:7.3f}".format, arr[i]))
            str_arr_i = ','.join(map(fmt.format, arr[i]))
            f.write(lst_ptn_names[i] + ";" + str_arr_i + "\n")

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
    for a,csta,cstb,ccta,cctb in zip(arr_list[slc_nwb],
                                     chunk_st_a[slc_nwb],
                                     chunk_st_b[slc_nwb],
                                     chunk_ct_a[slc_nwb],
                                     chunk_ct_b[slc_nwb]):
         mat_ret[csta:csta+ccta,cstb:cstb+cctb] = a.reshape(ccta,cctb)

    np.set_printoptions(edgeitems=30, linewidth=1000000, formatter=dict(float=lambda x: "%.3f" % x))
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


parser = argparse.ArgumentParser(description='rebuild matrix options')

parser.add_argument('-f', action="store", dest="sample_result_file")
parser.add_argument('-validate', action="store_true", default=False)
parser.add_argument('-csv', action="store_true", default=False)

args = parser.parse_args()

print("input_file = ",args.sample_result_file)
print("save csv = ",args.csv)
print("validate = ",args.validate)

f = h5py.File(args.sample_result_file + '.res_all.h5', 'r')
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

pn_filename = args.sample_result_file
if pn_filename.endswith('.h5'):
    pn_filename = pn_filename[:-3] + '.json'
print("protein_file_name:",pn_filename)
with open(pn_filename) as json_file:  
    lst_protein_names = json.load(json_file)
#print(lst_protein_names)

fmt_str="%7.3f" # "{:7.3f}"
fmt_str3="{:7.3f}"
if (args.csv):
    out_dir = args.sample_result_file + "_res_csv"
    mkdir_p(out_dir)
    print("Saving csv in:",out_dir)
    #np.savetxt(out_dir+"/normal.csv", mat_normal_h5, delimiter=",",fmt=fmt_str)
    save_mat(mat_normal_h5, fmt = fmt_str3, lst_ptn_names= lst_protein_names, out_csv_fname = out_dir+"/normal.csv")
    save_mat(mat_sarika_h5, fmt = fmt_str3, lst_ptn_names= lst_protein_names, out_csv_fname = out_dir+"/sarika.csv")
    save_mat(mat_generalised_h5, fmt = fmt_str3, lst_ptn_names= lst_protein_names, out_csv_fname = out_dir+"/generalised.csv")
    save_mat(mat_wu_h5, fmt = fmt_str3, lst_ptn_names= lst_protein_names, out_csv_fname = out_dir+"/wu.csv")
    #np.savetxt(out_dir+"/normal.csv", mat_normal_h5, delimiter=",",fmt=fmt_str)
    #np.savetxt(out_dir+"/sarika.csv", mat_sarika_h5, delimiter=",",fmt=fmt_str)
    #np.savetxt(out_dir+"/generalised.csv", mat_generalised_h5, delimiter=",",fmt=fmt_str)
    #np.savetxt(out_dir+"/wu.csv", mat_wu_h5, delimiter=",",fmt=fmt_str)

if (args.validate):
    sample_name = args.sample_result_file.split(".")[0]
    py_path_name = sample_name +"_csv_uint16/"
    print("Verifying results from:")
    print("py_path_name=",py_path_name)
    # for key in keys:
    mat_normal_py = np.loadtxt(py_path_name + "normal.csv",delimiter=",")
    mat_tol = 0.002
    ret = np.allclose(mat_normal_h5, mat_normal_py, atol=mat_tol)
    if (ret):
        print("normal is OK!")
    else:
        print("normal is not ok!")
        np.savetxt("normal.csv", mat_normal_h5, delimiter=",",fmt="%7.3f")
    
    # for key in keys:
    mat_sarika_py = np.loadtxt(py_path_name + "sarika.csv",delimiter=",")
    ret = np.allclose(mat_sarika_h5, mat_sarika_py, atol=mat_tol)
    if (ret):
        print("sarika is OK!")
    else:
        print("sarika is not ok!")
        #np.savetxt('sarika.csv',fmt="%7.3f")
        np.savetxt("sarika.csv", mat_sarika_h5, delimiter=",",fmt="%7.3f")
    
    mat_generalised_py = np.loadtxt(py_path_name + "generalised.csv",delimiter=",")
    ret = np.allclose(mat_generalised_h5, mat_generalised_py, atol=mat_tol)
    if (ret):
        print("generalised is OK!")
    else:
        print("generalised is not ok!")
        np.savetxt("generalised.csv", mat_generalised_h5, delimiter=",",fmt="%7.3f")
    
    mat_wu_py = np.loadtxt(py_path_name + "wu.csv",delimiter=",")
    ret = np.allclose(mat_wu_h5, mat_wu_py, atol=mat_tol)
    if (ret):
        print("wu is OK!")
    else:
        print("wu is not ok!")
        np.savetxt("wu.csv", mat_wu_h5, delimiter=",",fmt="%7.3f")

    #mat_cosine_py = np.loadtxt(py_path_name + "cosine.csv",delimiter=",")
    #ret = np.allclose(mat_cosine_h5, mat_cosine_py, atol=mat_tol)
    #if (ret):
    #    print("cosine is ok!")
    #else:
    #    print("cosine is not ok!")
    #print(mat_normal_py)
    #np.set_printoptions(edgeitems=30, linewidth=100000, formatter=dict(float=lambda x: "%7.3f" % x))
    #print(mat_cosine)
