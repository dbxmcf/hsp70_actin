
This README.md provides a brief introduction on how to use the distributed version to run 3D protein structure comparison in parallel. The code uses two modes of parallism by using the distributed memory programming model and shared memory model, namely:
- Massage Passing Interface (MPI) + Open Multi-Processing (OpenMP)
  - This version runs on regular multi-core nodes (without GPU)
- Massage Passing Interface (MPI) + Open Accelerators (OpenACC)
  - This version runs on GPU nodes 

The implementation strategy and background will be detailed in our later paper.

## Prerequisites/Software environment

The distributed version requires the following software environment to be loaded, python3 and pgi based mvapich2 and hdf5

```
module load python/3.5.2-anaconda-tensorflow
module load mvapich2.3/pgi-18.7
module load hdf5/mvapich2.3-pgi-18.7
```

## Clone the code repository from bitbucket:

From QB2 terminal, clone the code repository (with your bitbucket credential) using the below command:

```
    [<user>@qb2 <user>]$ pwd
    /work/<user>
    [<user>@qb4 <user>]$ git clone https://bitbucket.org/dbxmcf/hsp70_actin.git
    <use your bitbucket username/password>
```
Enter the repository folder:

```
    [<user>@qb4 <user>]$ cd hsp70_actin/
    [<user>@qb4 hsp70_actin]$
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

The protein structure file is large and consumes a lot of memory, according to current observations, a single line vector in the csv file (localFeatureVect_theta29_dist35_NoFeatureSelection_keyCombine0.csv) will take about 5MB memory storage (denote as `l`), therefore, a 1000 line sample will take about 5x1000=5GB in memory. Suppose we have an initial csv file of `N` lines, our current code will first divide the entire `N` lines into `C` chunks, each chunk will then contain `N/C` lines, each MPI process will take about `2xN/Cxl` memory, on QB2 workq compute nodes, each node has 64GB memory, the maximum number of csv lines that can be fit by a compute node should be less than ~64GB/(2x5MB) = 6000, in other words, we cannot use more than 6000 lines on a single compute node. To achieve the best performance, we suggest each MPI process does not exceed half of the available compute node memory, so each process should not use more than 3000 csv lines, i.e.:

max_csv_line_per_mpiprocess = 3000

### Number of MPI processes

In order to achieve a good load balancing, the number of line chunks `C` can only be an even number greater than 2 (using two chunks is computationally meaningless), so the possible values for `C` can only be 

4, 6, 8, 10,..., 

Correspondingly, the number of MPI processes \<NPROCS> we can use is based on the equation:

np = CxC/2

so the number of MPI processes for the above C values are:

8, 18, 32, 50,...,

### Number of compute nodes needed for a particular job

To calculate the number of processes and compute nodes that we should use for each particular problem, use the following rules:

1) Based on previous 3 example cases, in order to achieve the best performance, we need to use 2 MPI processes per compute node, i.e., use "-ppn 2" in the mpirun command.
2) Based on 1), we can calculate the number of compute nodes needed to use for each value of C is:
    num_nodes = CxC/2/ppn = CxC/2/2 = CxC/4
3) For above C values, the numbers of compute nodes needed are: 4, 9, 16, 25,...,

### Maximum size of sample we can calculate for each type of job

Using the above setting, the maximum lines of csv file a job can process can be calculated as:

total_csv_lines = max_csv_line_per_mpiprocess x ppn x num_nodes

As an example, for num_nodes = 4, 

total_csv_lines = 3000 x 2 x 4 = 24000




### 

### How to prepare input file

The input file is generated by running the lf_csv2_dtype_h5.py from the terminal.

```
$ python lf_csv2_dtype_h5.py -f sample_hsp70_actin [-o out_dir]
```

The above python script will convert the last step csv file () converted to an hdf5 file named `"out_dir/sample_hsp70_actin.h5"` to be used as the input file for the distributed version of the program.

### How to run the program

To run the program, use the mpirun command with the following syntax:

```
$ mpirun -np <NPROCS> -ppn <processes_per_node> ./pomp -f <inputfile_name>
```

So in order to run the analysis using sample_hsp70_actin with 8 processes, the correct command to use is:

```
$ mpirun -np 8 -ppn 2 ./pomp -f sample_hsp70_actin.h5
```

## Distribute the workload to different compute nodes: How to combine the \<NPROCS\> and -ppn values: 

To evenly distribute the workload to different compute nodes, use the -ppn <processes_per_node>

### How to process the output file

As mentioned earlier, the output of the program is also an hdf5 file with all the result information from each MPI process, however the calculated results from each MPI process are not stored according to the original order of the matrix, a python code will be used to re-organize and gather information from the resultant hdf5 file and combine them into a numpy readable matrix (2D table) and can be used for subsequent analysis. To run the output-processing script, use the below syntax:

```
python rebuild_mat.py -f $input_h5 -csv
```





