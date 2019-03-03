import numpy as np
size = 20000
a = np.random.random_sample((size, size),dtype=bool)
b = np.random.random_sample((size, size),dtype=bool)
#n = np.dot(a,b)
n = a & b
