#!/usr/bin/env python
from mpi4py import MPI
import numpy as np

# https://stackoverflow.com/questions/36025188/along-what-axis-does-mpi4py-scatterv-function-split-a-numpy-array/36082684#36082684
# https://stackoverflow.com/questions/6081008/dump-a-numpy-array-into-a-csv-file?rq=1
# https://info.gwdg.de/~ceulig/docs-dev/doku.php?id=en:services:application_services:high_performance_computing:mpi4py
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

sendbuf = np.zeros(10, dtype='i') + rank
recvbuf = None
if rank == 0:
    recvbuf = np.empty([size, 10], dtype='i')
comm.Gather(sendbuf, recvbuf, root=0)
if rank == 0:
    for i in range(size):
        #assert np.allclose(recvbuf[i,:], i)
        print(recvbuf[i,:], i)

