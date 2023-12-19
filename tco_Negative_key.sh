#!/bin/bash
#PBS -l nodes=4:ppn=20
#PBS -l walltime=01:00:00
#PBS -q checkpt
#PBS -N jcd_mpi_omp
#PBS -j oe
#-A hpc_hpcadmin5
#PBS -A loni_loniadmin1

#cd $PBS_O_WORKDIR
#module purge
#module load python/3.6.2-anaconda-tensorflow
#module load hdf5/mvapich2.3-pgi-18.7

# fill in absolute path here for input, *must* contain the csv file
module purge
module load python/3.7.6
module load mvapich2/2.3.3/intel-19.0.5-hydra
cd ../
SECONDS=0
python aa_Triplets_key_cal_Negative_key.py
echo "took $SECONDS sec."
cd -
#input_sample_folder="./sample_hsp70_actin"
input_sample_folder="../Negative_key"
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

#export OMP_NUM_THREADS=1
#SECONDS=0
#mpirun -np $NPROCS -f $PBS_NODEFILE -ppn 2 ./pgi.mpi.pacc.out -f $input_h5
#echo "gpu run $i GPU, $SECONDS sec"

python rebuild_mat.py -f $input_h5 -csv #-validate

