
This README.md provides a brief introduction on how to use the distributed version to run 3D protein structure comparison in parallel. The code uses two modes of parallism by using the distributed memory programming model and shared memory model, namely:
- Massage Passing Interface (MPI) + Open Multi-Processing (OpenMP)
  - This version runs on regular multi-core nodes (without GPU)
- Massage Passing Interface (MPI) + Open Accelerators (OpenACC)
  - This version runs on GPU nodes 

The implementation strategy and background will be detailed in the XXX paper, 

## To the impatient

Jump to the section <a href="Example script">Example script</a>, which should allow running a job without going through all the details

## How to compile the executable file

1. CPU Version

To compile the MPI+OpenMP executable, type the below command:

 if only one version is needed, use the command:

```sh
$ make pomp
```

2. GPU Version

The GPU version is typically about twice as fast as the CPU multi-threaded version.

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

### 

### How to prepare input file

The input file is generated by running the lf_csv2_dtype_h5.py from the terminal.

```bash
$ python lf_csv2_dtype_h5.py -f sample_hsp70_actin [-o out_dir]
```

The above python script will convert the last step csv file () converted to an hdf5 file named `"out_dir/sample_hsp70_actin.h5"` to be used as the input file for the distributed version of the program.

### How to run the program

To run the program, use the mpirun command with the following syntax:

```bash
$ mpirun -np <NPROCS> -ppn <processes_per_node> ./pomp -f <inputfile_name>
```

So in order to run the analysis using sample_hsp70_actin with 8 processes, the correct command to use is:

```sh
$ mpirun -np 8 -ppn 2 ./pomp -f sample_hsp70_actin.h5
```

## Distribute the workload to different compute nodes: How to combine the \<NPROCS\> and -ppn values: 

To evenly distribute the workload to different compute nodes, use the -ppn <processes_per_node>

### How to process the output file

As mentioned earlier, the output of the program is also an hdf5 file with all the result information from each MPI process, however the calculated results from each MPI process are not stored according to the original order of the matrix, a python code will be used to re-organize and gather information from the resultant hdf5 file and combine them into a numpy readable matrix (2D table) and can be used for subsequent analysis. To run the output-processing script, use the below syntax:

```bash
python rebuild_mat.py -f $input_h5 -csv
```

## Example script

1. Copy the template job script to another name:

    ```bash
    cp template_large.pbs my_script.pbs
    ```

2. Change line 15, input_sample_folder to the desired folder name, this sample folder *must* contain the input csv file `"<sample_folder_root>/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv"` from the last step, e.g.:

    ```
    15 input_sample_folder="sample_protease_mix_1"
    ```
3. Change line 16, output_folder to a desired output folder name, if the output folder does not exist, the folder will be created, e.g.:
    ```
    16 output_folder="out"
    ```
4. Change the allocation name to your allocation:

    ```
    #PBS -A <your_allocation_name>
    ```



