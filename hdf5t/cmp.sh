#!/bin/bash
make pgiacc
sample_name="sample_a-b_mix_2"
mpirun -np 8 ./pgiacc.mpi.acc.out -f ${sample_name}.h5
python ./rebuild_mat.py ${sample_name}

