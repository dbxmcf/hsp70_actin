
# Introduction

This README.md provides a brief introduction on how to use the distributed version to run 3D protein structure comparison in parallel. The code uses two modes of parallism by using the distributed memory programming model and shared memory model, namely:

- Massage Passing Interface (MPI) + Open Multi-Processing (OpenMP)
    - This version runs on regular multi-core nodes (without GPU)

- Massage Passing Interface (MPI) + Open Accelerators (OpenACC)
    - This version runs on GPU nodes 

The implementation strategy and background will be detailed in our later paper.

## Prerequisites/Software environment

The distributed version requires the following software environment to be loaded, python3 and pgi based mvapich2 and hdf5, to avoid conflict between other modules, do a `module purge` before loading the following modules:

```
module purge
module load python/3.5.2-anaconda-tensorflow
module load mvapich2.3/pgi-18.7
module load hdf5/mvapich2.3-pgi-18.7
```

## Clone the code repository from bitbucket:

From QB2 terminal, clone the code repository (with your bitbucket credential) using the below command:

```
    [<user>@qb2 <user>]$ pwd
    /work/<user>
    [<user>@qb2 <user>]$ git clone https://bitbucket.org/dbxmcf/hsp70_actin.git
    <use your bitbucket username/password>
```
Enter the repository folder:

```
    [<user>@qb2 <user>]$ cd hsp70_actin/
    [<user>@qb2 hsp70_actin]$
```

## Compile the executable file

1. CPU Version

To compile the MPI+OpenMP executable, type the below command:

 if only one version is needed, use the command:

```sh
$ make pomp
```

1. GPU Version

The GPU version is typically about twice as fast as the CPU multi-threaded version however subject to memory limitation at this moment:

To compile the MPI+OpenACC executable, type the below command:

```sh
$ make pacc
```

1. Compile both versions

To compile both versions (CUP and GPU) using one command, type:

```sh
$ make all
```

CPU version should be sufficient for a typical job less than 50k lines of csv input.

## Example job scripts

1. Copy the template job script to another name, e.g.:

    ```
    cp template_large.pbs my_script.pbs
    ```

1. Open the new script and change line 15, input_sample_folder to the desired input folder name, note that this sample folder *must* contain the input csv file 
``<sample_folder_root>/theta29_dist35/localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv`` 
from the last step, e.g.:

    ```
    input_sample_folder="sample_protease_mix_1"
    ```

1. Change line 16, output_folder to a desired output folder name, if the output folder does not exist, the folder will be created, e.g.:

    ```
    output_folder="/work/<user>/sim_results"
    ```

1. Also remember to change the allocation name to your allocation:

    ```
    PBS -A <your_allocation_name>
    ```

1. Submit the job

    ```
    qsub my_script.pbs
    ```

### Memory usage consideration

The 3D protein structure file is large and consumes a lot of memory, according to current observations, a single line vector in the csv file (localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv) will take about 5MB memory storage (denote as `l`), therefore, a 1000 line sample will take about 5x1000=5GB in memory. Suppose we have an initial csv file of `N` lines, our current code will first divide the entire `N` lines into `C` chunks, each chunk will then contain `N/C` lines, each MPI process will take two chunks of data for comparison which takes about `2xN/Cxl` memory, on QB2 workq compute nodes, each node has 64GB memory, the maximum number of csv lines that can fit within a compute node should be less than ~64GB/(2x5MB) = 6400. In other words, we cannot use more than 6000 lines on a single compute node. Further, to achieve the best performance, we suggest each MPI process does not exceed half of the available compute node memory (so that one compute node can fit two MPI processes) which means each process does not exceed ~3200 csv lines, i.e.:

max_csv_line_per_mpiprocess = 3200

### Number of MPI processes

In order to achieve a good load balancing, the number of line chunks `C` is currently limited to even numbers (using two chunks is still acceptable but might result in error in some cases), so the possible values for `C` can only be 

4, 6, 8, 10,..., 

Correspondingly, the number of MPI processes \<NPROCS> we can use is based on below equation:

np = CxC/2 (details on derivation will be in our paper)

so the number of MPI processes for the above C values are:

8, 18, 32, 50,...,

### Number of compute nodes needed for a particular job

To calculate the number of processes and compute nodes that we should use for each job, use the following rules:

1) Based on previous 3 example cases, in order to achieve the best performance, we need to use 2 MPI processes per compute node, i.e., use "-ppn 2" in the *mpirun* command. (Note: Please do *not* confuse this option with the ppn on the top part of the PBS job script, in nodes=\<N>:ppn=20 means 20 cores per node, they are the same name however completely different meanings!)
2) Based on 1), we can calculate the number of compute nodes needed to use for each value of C is:
    num_nodes = CxC/2/ppn = CxC/2/2 = CxC/4
3) For above C values, the numbers of compute nodes needed are: 4, 9, 16, 25,...,

### Maximum size of sample (total number of lines in the csv file) we can process for each job

Using the above setting, the maximum lines of csv file a job can process can be calculated as:

total_csv_lines = max_csv_line_per_mpiprocess x C

As used in the example job script, for C = 4, np = 8, num_nodes = 4,

total_csv_lines = 3200 x 4 = 12800

It follows that, if we ask for 9 nodes in our job script:


```
PBS nodes=9:ppn=20
```


With -ppn=2 in the mpirun command, this corresponds to C = sqrt(2 x np x ppn) = sqrt(2x9x2) = sqrt(36) = 6, and:

total_csv_lines = 3200 x 6 = 19200

For C=32, np = 512, num_nodes = 256 (maximum allowed nodes for a single job):

total_csv_lines = 3200 x 32 = 102400

Thus to *get the best performance*, 100k lines of csv is about the limit of a single job on QB2. 

Below table summarized the number of data chunks `C`， the needed MPI processes for `C`， number of compute nodes needed based on ppn=2 and corresponding maximum csv lines a job is able to process:


| C  | np (number of MPI processes) | num_nodes (based on -ppn=2) | max csv lines |
|----|------------------------------|----------------------------|---------------|
| 2  | 2                            | 1                          | 6400          |
| 4  | 8                            | 4                          | 12800         |
| 6  | 18                           | 9                          | 19200         |
| 8  | 32                           | 16                         | 25600         |
| 10 | 50                           | 25                         | 32000         |
| 12 | 72                           | 36                         | 38400         |
| 14 | 98                           | 49                         | 44800         |
| 16 | 128                          | 64                         | 51200         |
| 18 | 162                          | 81                         | 57600         |
| 20 | 200                          | 100                        | 64000         |
| 22 | 242                          | 121                        | 70400         |
| 24 | 288                          | 144                        | 76800         |
| 26 | 338                          | 169                        | 83200         |
| 28 | 392                          | 196                        | 89600         |
| 30 | 450                          | 225                        | 96000         |
| 32 | 512                          | 256                        | 102400        |

**Note** The above table is based on two MPI processes per node (-ppn=2), if we use only one process on each compute node, i.e., ppn=1, then each node can process 6400 lines, the maximum allowable csv lines is approximately:

total_csv_lines = 6400 x 22 = 140800

### How to prepare input file (You can ignore the below contents at this moment if you are using the template job script)

The input file is generated by running the lf_csv2_dtype_h5.py from the terminal.

```
$ python lf_csv2_dtype_h5.py -f sample_hsp70_actin [-o out_dir]
```

The above python script will convert the last step csv file converted to an hdf5 file named `"out_dir/sample_hsp70_actin.h5"` to be used as the input file for the distributed version of the program.

### How to run the program

To run the program, use the mpirun command with the following syntax:

```
$ mpirun -np <NPROCS> -ppn <processes_per_node> ./pomp -f <inputfile_name>
```

So in order to run the analysis using sample_hsp70_actin with 8 processes, the correct command to use is:

```
$ mpirun -np 8 -ppn 2 ./pomp -f sample_hsp70_actin.h5
```

### Distribute the workload to different compute nodes: How to combine the \<NPROCS\> and -ppn values: 

To evenly distribute the workload to different compute nodes, use the -ppn <processes_per_node>

### How to process the output file

As mentioned earlier, the output of the program is also an hdf5 file with all the result information from each MPI process, however the calculated results from each MPI process are not stored according to the original order of the matrix, a python code will be used to re-organize and gather information from the resultant hdf5 file and combine them into a numpy readable matrix (2D table) and can be used for subsequent analysis. To run the output-processing script, use the below syntax:

```
python rebuild_mat.py -f $input_h5 -csv
```





