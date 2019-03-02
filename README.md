
This README.md provides a brief introduction on how to use the distributed version to run 3D protein structure comparison in parallel. The code uses two modes of parallism by using the distributed memory programming model and shared memory model, namely:
- Massage Passing Interface (MPI) + Open Multi-Processing (OpenMP)
  - This version runs on regular multi-core nodes (without GPU)
- Massage Passing Interface (MPI) + Open Accelerators (OpenACC)
  - This version runs on GPU nodes 

The implementation strategy and background will be detailed in the XXX paper, 

## How to compile the executable file

1. CPU Version

To compile the MPI+OpenMP executable, type the below command:

 if only one version is needed, use the command:

```sh
$ make pomp
```

2. GPU Version

To compile the MPI+OpenACC executable, type the below command:

```sh
$ make pacc
```

This will compile all two versions of the code, the MPI+OpenMP version can run on multi-core node (without GPU), 

to compile the MPI + OpenMP version

Massage Passing Interface (MPI) + Open Multi-Processing (OpenMP)

## How to run jobs

### Prerequisites

#### Software environment

### Memory usage consideration

The protein structure file is large and consumes a lot of memory when loaded to par_hybrid_cmp, according to current observations, a 

### Usage of hdf5

Using hdf5 will greatly reduce the file size and also the speed when doing parallel input/output, therefore we use the hdf5 (.h5) as the input file for the protein structure and also for the output matrices (normal, generalises, wu and sarika).

### How to prepare input file

The input file is generated by running the below command in the terminal.

```python
$ python lf_csv2hdf5.py sample_hsp70_actin
```

The previous csv file will be converted to an hdf5 file named sample_hsp70_actin.h5 and will be used as the input file for the program.

### How to run the program

To run the program, use the mpirun command with the following syntax:

```sh
$ mpirun -np <NPROCS> ./pomp -f <inputfile_name>
```

So in order to run the analysis using sample_hsp70_actin with 8 processes, the correct command to use is:

```sh
$ mpirun -np 8 ./pomp -f sample_hsp70_actin.h5
```

### Distribute the workload to different compute nodes

The primary feature of the code is its capability of distributing the workload (in this case, high memory usage) to different compute nodes so that large size comparison can be made possible. For a fixed given sample size, the more number of processes used, the less memory will be allocated to each process. For a fixed size memory compute node, when more compute nodes are used, the larger the total sample size can be analyzed. 

### How to use output file

As mentioned earlier, the output of the program is also an hdf5 file with all the result information from each MPI process, however the calculated results from each MPI process are not stored according to the original order of the matrix, a python code will be used to gather information from the result hdf5 file and combine them into a numpy readable matrix (2D table) and can be used for subsequent analysis. 


