#!/usr/bin/env python
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

a_size = 4
recvdata = None
senddata = (rank+1)*np.arange(a_size,dtype=np.float64)
if rank == 0:
   recvdata = np.arange(size*a_size,dtype=np.float64)
comm.Gather(senddata,recvdata,root=0)
print('on task',rank,'after Gather:    data = ',recvdata)

counts=(2,3,4)
dspls=(0,3,8)
if rank == 0:
   recvdata = np.empty(12,dtype=np.float64)
sendbuf = [senddata,counts[rank]]
recvbuf = [recvdata,counts,dspls,MPI.DOUBLE]
comm.Gatherv(sendbuf,recvbuf,root=0)
print('on task',rank,'after Gatherv:    data = ',recvdata)
