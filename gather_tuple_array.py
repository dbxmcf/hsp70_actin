#!/usr/bin/env python
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

m_idx_type = np.int
mpi_type = MPI.LONG

a_size = 6
recvdata = None

#list_rank = 

#senddata = (rank+1)*np.arange(a_size,dtype=m_idx_type)
#senddata = (rank+1)*np.arange(a_size,dtype=m_idx_type)
senddata = (rank+1)*np.arange(a_size,dtype=m_idx_type) #.reshape(-1,2)

if rank == 2:
    senddata = (rank+1)*np.arange(a_size/2+1,dtype=m_idx_type)
print('On rank', rank, 'senddata=', senddata)

#if rank == 0:
#   recvdata = np.arange(size*a_size,dtype=m_idx_type)
#comm.Gather(senddata,recvdata,root=0)
#print('on task',rank,'after Gather:    data = ',recvdata)

#counts=(2,3,4)
counts = np.zeros(size)
dspls = np.zeros_like(counts)
counts[rank] = len(senddata)

rank_div = 2

if rank <=rank_div:
   dspls[rank] = rank*len(senddata)
else:
   dspls[rank] = rank*
dspls[rank] = (0,3,8)
if rank == 0:
   #recvdata = np.empty(12,dtype=m_idx_type)
   recvdata = np.zeros(12,dtype=m_idx_type)
sendbuf = [senddata,counts[rank]]
recvbuf = [recvdata,counts,dspls,mpi_type]
comm.Gatherv(sendbuf,recvbuf,root=0)
if rank ==0:
   print('on task',rank,'after Gatherv:    data = ',recvdata)
