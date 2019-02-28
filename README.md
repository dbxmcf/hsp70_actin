# Instruction

This README.md provides a brief introduction on how to use the distributed version to run 3D protein structure in parallel. The code uses two modes of hybrid - parallism, namely:
- Massage Passing Interface (MPI) + Open Multi-Processing (OpenMP)
- This version runs- 
- Massage Passing Interface (MPI) + Open Accelerators (OpenACC)

## How to compile the executable file

To compile the executable, type the below command:

```sh
$ make all
```
This will compile all two versions of the code, the MPI+OpenMP version can run on multi-core node (without GPU), 

 if only one version is needed, use the command:

```sh
$ make pomp
```

to compile the MPI + OpenMP version

Massage Passing Interface (MPI) + Open Multi-Processing (OpenMP)

## How to run jobs

### Memory usage consideration

The protein structure file is large and consumes a lot of memory when loaded to the program

### How to prepare input file

### How to use output file

#### Usage of hdf5

Using hdf5 will greatly reduce the file size and also the speed when doing parallel input/output, therefore we use the hdf5 (.h5) as the input file for the protein structure and also for the output matrices (normal, generalises, wu and sarika).





