#!/usr/bin/env python

from __future__ import print_function
import sys
import errno 
import os
import time
import numpy as np
import argparse
import json
#from numpy import linalg as LA
#import StringIO
#import itertools
#from scipy import spatial
#import cartesian
#from cartesian import *
#from itertools import combinations
import h5py

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

parser = argparse.ArgumentParser(description='Generate hdf5 from csv file:')
#parser.add_argument('-f', action="store", dest="sample_name", required=True, help='Input folder name')
parser.add_argument('-f','--input_folder', action="store", required=True, help='Input folder name')
parser.add_argument('-o','--output_folder', action="store",required=False, 
                    help='Output folder name', default='hdf5t')
#args = parser.parse_args(['-h'])
args = parser.parse_args()


#fname = "test.csv"
#fname = "sample_hsp70_actin/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"
#fname = sample_name + "/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"
#fname = "sample_protease_mix_1/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"

#sample_name = "sample_hsp70_actin"
#sample_name = "sample_a-b_mix_2"
#sample_name = "sample_protease_mix_1"
sample_folder = args.input_folder
sample_folder = sample_folder.rstrip('//')
#get the last path as the sample name
sample_name = os.path.basename(os.path.normpath(sample_folder))
print("input_folder = ",sample_folder)
fname = sample_folder + "/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"

#sample_name = "hdf5t"
#fname = sample_name + "/ta.csv"

start_time=time.time()

arrs = []
m_datatype = np.uint16
protein_names = []
#m_datatype = np.int64

with open(fname) as fcsv:
    lines=fcsv.readlines()
    #for idx,line in enumerate(lines[n_st:]):
    for idx,line in enumerate(lines):
        l = list(line.split(';')[1].split(','))
        #l_arr = np.asarray(l[:-1]).astype(np.float) 
        l_arr = np.asarray(l[:-1],dtype=m_datatype)
        #l_arr = np.asarray(l[:],dtype=m_datatype)
        arrs.append(l_arr)
        protein_names.append(list(line.split(';'))[0])

#print(protein_names)
data = np.array(arrs,dtype=m_datatype)
#print(data)
print(data.shape)
print(data.dtype)
print("min val=",data.min())
print("max val=",data.max())
num_non_zeros=np.count_nonzero(data)
print("non_zeros=",num_non_zeros)
print("non_zeros_percent=",num_non_zeros/data.shape[0]/data.shape[1])

end_time=time.time()
total_time=((end_time)-(start_time))
print("Time taken for reading csv: {}".format(total_time))

#h5_filename = "hdf5t/" + sample_name + "_dtype.h5"
mkdir_p(args.output_folder)
filename_prefix = args.output_folder + "/" + sample_name

h5_filename = filename_prefix + ".h5"
h5f = h5py.File(h5_filename, 'w')
h5f.create_dataset('Data1', data=data, dtype=m_datatype)
ascii_protein_names = [n.encode("ascii", "ignore") for n in protein_names]
h5f.create_dataset('ProteinNames', (len(ascii_protein_names),1),'S10', ascii_protein_names)
#h5f.create_dataset('ProteinNames', data=protein_names, dtype=h5py.special_dtype(vlen=str))
h5f.close()
print(h5_filename + " file created.")

pn_filename = filename_prefix + ".json"
with open(pn_filename, 'w') as outfile:
    json.dump(protein_names, outfile)
print(pn_filename + " file created.")

