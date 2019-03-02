import numpy as np
from itertools import combinations, chain
from scipy.special import comb

def comb_index(n, k):
    count = comb(n, k, exact=True)
    index = np.fromiter(chain.from_iterable(combinations(range(n), k)), 
            int, count=count*k)
    return index.reshape(-1, k)

#data = np.array([[1,2,3,4,5],[10,11,12,13,14]])
data = np.array([[1,2,3],[1,2,3]])
idx = comb_index(3, 2)
print(data[:, idx])
