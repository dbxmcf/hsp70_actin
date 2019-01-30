#!/usr/bin/env python

#!/usr/bin/env python
import numpy as np
import pandas as pd
import h5py

val_range = 16
n = 17
m = 12
#test_arr = np.random.randint(val_range, size=(n, m),dtype=np.uint16)
data_type = np.int
test_arr = np.zeros((n,m),dtype=data_type)
test_arr[0:n,:]=np.arange(n).reshape(n,-1)*100
test_arr += np.arange(m,dtype=data_type)
h5f = h5py.File('ta.h5', 'w')
h5f.create_dataset('Data1', data=test_arr)
h5f.close()

print(test_arr)
