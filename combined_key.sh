#!/bin/bash

# fill in absolute path here for input, *must* contain the csv file
#module purge
#module load python/3.7.6
#module load mvapich2/2.3.3/intel-19.0.5-hydra
#cd ../
# build the executable
make -f makefile.mpich
SECONDS=0
python aa_Triplets_key_cal_Negative_key.py
echo "took $SECONDS sec."
#cd -
#input_sample_folder="./sample_hsp70_actin"
input_sample_folder="./Negative_key"
# fill in output path here
output_folder="./hdf5t"

# remove possible tailing slashes
output_folder=${output_folder%/}
mkdir -p ${output_folder}
input_sample_folder=${input_sample_folder%/}

python lf_csv2_dtype_h5.py -f "$input_sample_folder" -o "$output_folder"

input_sample_name=${input_sample_folder##*/}
echo "input_sample_name=$input_sample_name"

input_h5="${output_folder}/${input_sample_name}.h5"
echo "input_h5=$input_h5"

echo "using $NPROCS processes..."

export OMP_NUM_THREADS=1
SECONDS=0
#mpirun -np $NPROCS -f $PBS_NODEFILE -ppn 2 ./pgi.mpi.pomp.out -f $input_h5
export SLURM_OVERLAP=1
#mpirun -np 8 ./mpichomp.out -f $input_h5
mpirun -np 8 ./mpichomp.out -f $input_h5
echo "omp run trd=$OMP_NUM_THREADS, $SECONDS sec"

python rebuild_mat.py -f $input_h5 -csv #-validate

# fill in absolute path here for input, *must* contain the csv file
SECONDS=0
python aa_Triplets_key_cal_without_Negative_key.py
echo "took $SECONDS sec."

input_sample_folder="./without_Negative_key/"
# fill in absolute path here for output
output_folder="./hdf5t"

# remove possible tailing slashes
output_folder=${output_folder%/}
input_sample_folder=${input_sample_folder%/}

python lf_csv2_dtype_h5.py -f "$input_sample_folder" -o "$output_folder"

input_sample_name=${input_sample_folder##*/}
echo "input_sample_name=$input_sample_name"

input_h5="${output_folder}/${input_sample_name}.h5"
echo "input_h5=$input_h5"

echo "using $NPROCS processes..."

export OMP_NUM_THREADS=1
SECONDS=0
#mpirun -np $NPROCS -f $PBS_NODEFILE -ppn 2 ./pgi.mpi.pomp.out -f $input_h5
export SLURM_OVERLAP=1
#mpirun -np 8 ./mpichomp.out -f $input_h5
mpirun -np 8 ./mpichomp.out -f $input_h5
echo "omp run trd=$OMP_NUM_THREADS, $SECONDS sec"

python rebuild_mat.py -f $input_h5 -csv #-validate

