#!/usr/bin/env python
import numpy as np
#import pandas as pd
import h5py

val_range = 1024
#n = 4295
#n = 16000
#n = 70000
#n = 100000
#n = 160000
n = 20000
#m = 1417419
m = 1500000
#m = 1000000
m_type = np.uint16
test_arr = np.random.randint(val_range, size=(n, m),dtype=m_type)

#data_type = np.int
#test_arr = np.zeros((n,m),dtype=data_type)
#test_arr[0:n,:]=np.arange(n).reshape(n,-1)*100
#test_arr += np.arange(m,dtype=data_type)
#np.savetxt("ta.csv", test_arr, delimiter=",",fmt='%5d')

fname = 'asdf/ta'+str(n)+'.h5'
h5f = h5py.File(fname, 'w')
h5f.create_dataset('Data1', data=test_arr,dtype=m_type)
h5f.close()

#print(test_arr)
