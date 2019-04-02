/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the COPYING file, which can be found at the root of the source code       *
 * distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
 * If you do not have access to either file, you may request a copy from     *
 * help@hdfgroup.org.                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Example of using the parallel HDF5 library to access datasets.
 * Last revised: April 24, 2001.
 *
 * This program contains two parts.  In the first part, the mpi processes
 * collectively create a new parallel HDF5 file and create two fixed
 * dimension datasets in it.  Then each process writes a hyperslab into
 * each dataset in an independent mode.  All processes collectively
 * close the datasets and the file.
 * In the second part, the processes collectively open the created file
 * and the two datasets in it.  Then each process reads a hyperslab from
 * each dataset in an independent mode and prints them out.
 * All processes collectively close the datasets and the file.
 *
 * The need of requirement of parallel file prefix is that in general
 * the current working directory in which compiling is done, is not suitable
 * for parallel I/O and there is no standard pathname for parallel file
 * systems.  In some cases, the parallel file name may even needs some
 * parallel file type prefix such as: "pfs:/GF/...".  Therefore, this
 * example requires an explicite parallel file prefix.  See the usage
 * for more detail.
 */

#ifdef _OPENACC
#include <accel.h>              // OpenACC
#endif
#include <assert.h>
#include "hdf5.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "dynamic_2d_array.h"
#include "pair_cmbs.h"
#include "calc_coeffs.h"

#ifdef H5_HAVE_PARALLEL
/* Temporary source code */
#define FAIL -1
/* temporary code end */

/* Define some handy debugging shorthands, routines, ... */
/* debugging tools */
#define MESG(x)\
    if (verbose) printf("%s\n", x);\

#define MPI_BANNER(mesg)\
{printf("--------------------------------\n");\
    printf("Proc %d: ", mpi_rank); \
    printf("*** %s\n", mesg);\
    printf("--------------------------------\n");}

#define SYNC(comm)\
{MPI_BANNER("doing a SYNC"); MPI_Barrier(comm); MPI_BANNER("SYNC DONE");}
/* End of Define some handy debugging shorthands, routines, ... */

/* Constants definitions */
/* 24 is a multiple of 2, 3, 4, 6, 8, 12.  Neat for parallel tests. */
#define SPACE1_DIM1	24
#define SPACE1_DIM2	24
#define SPACE1_RANK	2
#define DATASETNAME1	"Data1"
#define DATASETNAME2	"Data2"
#define DATASETNAME3	"Data3"
/* hyperslab layout styles */
#define BYROW		1	/* divide into slabs of rows */
#define BYCOL		2	/* divide into blocks of columns */

#define PARAPREFIX	"HDF5_PARAPREFIX"	/* file prefix environment variable name */


/* dataset data type.  Int's can be easily octo dumped. */
typedef int DATATYPE;
//typedef short DATATYPE;

/* global variables */
int nerrors = 0;				/* errors count */
#ifndef PATH_MAX
#define PATH_MAX    512
#endif  /* !PATH_MAX */
char    testfiles[2][PATH_MAX];


int mpi_size, mpi_rank;				/* mpi variables */
int gpunum=0,ngpus=0;

/* option flags */
int verbose = 0;			/* verbose, default as no. */
int debug_mpi_rank = 0;		/* specify an mpi rank to print */
int debug_info = 0;          /* enable print mpi infos*/
int doread=1;				/* read test */
int dowrite=0;				/* write test */
int docleanup=0;			/* cleanup */

/* Prototypes */
void slab_set(hsize_t start[], hsize_t count[], hsize_t stride[], int mode);
void dataset_fill(hsize_t start[], hsize_t count[], hsize_t stride[], DATATYPE * dataset);
void dataset_print(hsize_t start[], hsize_t count[], hsize_t stride[], DATATYPE * dataset);
int dataset_vrfy(hsize_t start[], hsize_t count[], hsize_t stride[], DATATYPE *dataset, DATATYPE *original);
void phdf5writeInd(char *filename);
void phdf5readInd(char *filename);
void phdf5writeAll(char *filename, result_pointers_diagnol* result_data);
void phdf5readAll(char *filename);
void test_split_comm_access(char filenames[][PATH_MAX]);
int parse_options(int argc, char **argv);
void usage(void);
int mkfilenames(char *prefix);
void cleanup(void);


/*
 * Setup the dimensions of the hyperslab.
 * Two modes--by rows or by columns.
 * Assume dimension rank is 2.
 */
    void
slab_set(hsize_t start[], hsize_t count[], hsize_t stride[], int mode)
{
    switch (mode){
        case BYROW:
            /* Each process takes a slabs of rows. */
            stride[0] = 1;
            stride[1] = 1;
            count[0] = SPACE1_DIM1/mpi_size;
            count[1] = SPACE1_DIM2;
            start[0] = mpi_rank*count[0];
            start[1] = 0;
            break;
        case BYCOL:
            /* Each process takes a block of columns. */
            stride[0] = 1;
            stride[1] = 1;
            count[0] = SPACE1_DIM1;
            count[1] = SPACE1_DIM2/mpi_size;
            start[0] = 0;
            start[1] = mpi_rank*count[1];
            break;
        default:
            /* Unknown mode.  Set it to cover the whole dataset. */
            printf("unknown slab_set mode (%d)\n", mode);
            stride[0] = 1;
            stride[1] = 1;
            count[0] = SPACE1_DIM1;
            count[1] = SPACE1_DIM2;
            start[0] = 0;
            start[1] = 0;
            break;
    }
}


/*
 * Fill the dataset with trivial data for testing.
 * Assume dimension rank is 2 and data is stored contiguous.
 */
    void
dataset_fill(hsize_t start[], hsize_t count[], hsize_t stride[], DATATYPE * dataset)
{
    DATATYPE *dataptr = dataset;
    hsize_t i, j;

    /* put some trivial data in the data_array */
    for (i=0; i < count[0]; i++){
        for (j=0; j < count[1]; j++){
            *dataptr++ = (i*stride[0]+start[0])*100 + (j*stride[1]+start[1]+1);
        }
    }
}


/*
 * Print the content of the dataset.
 */
void dataset_print(hsize_t start[], hsize_t count[], hsize_t stride[], DATATYPE * dataset)
{
    DATATYPE *dataptr = dataset;
    hsize_t i, j;

    /* print the slab read */
    for (i=0; i < count[0]; i++){
        printf("Row %lu: ", (unsigned long)(i*stride[0]+start[0]));
        for (j=0; j < count[1]; j++){
            printf("%03d ", *dataptr++);
        }
        printf("\n");
    }
}


/*
 * Print the content of the dataset.
 */
int dataset_vrfy(hsize_t start[], hsize_t count[], hsize_t stride[], DATATYPE *dataset, DATATYPE *original)
{
#define MAX_ERR_REPORT	10		/* Maximum number of errors reported */

    hsize_t i, j;
    int nerr;

    /* print it if verbose */
    if (verbose)
        dataset_print(start, count, stride, dataset);

    nerr = 0;
    for (i=0; i < count[0]; i++){
        for (j=0; j < count[1]; j++){
            if (*dataset++ != *original++){
                nerr++;
                if (nerr <= MAX_ERR_REPORT){
                    printf("Dataset Verify failed at [%lu][%lu](row %lu, col %lu): expect %d, got %d\n",
                            (unsigned long)i, (unsigned long)j,
                            (unsigned long)(i*stride[0]+start[0]), (unsigned long)(j*stride[1]+start[1]),
                            *(dataset-1), *(original-1));
                }
            }
        }
    }
    if (nerr > MAX_ERR_REPORT)
        printf("[more errors ...]\n");
    if (nerr)
        printf("%d errors found in dataset_vrfy\n", nerr);
    return(nerr);
}


/*
 * Example of using the parallel HDF5 library to create two datasets
 * in one HDF5 files with parallel MPIO access support.
 * The Datasets are of sizes (number-of-mpi-processes x DIM1) x DIM2.
 * Each process controls only a slab of size DIM1 x DIM2 within each
 * dataset.
 */

    void
phdf5writeInd(char *filename)
{
    hid_t fid1;			/* HDF5 file IDs */
    hid_t acc_tpl1;		/* File access templates */
    hid_t sid1;   		/* Dataspace ID */
    hid_t file_dataspace;	/* File dataspace ID */
    hid_t mem_dataspace;	/* memory dataspace ID */
    hid_t dataset1, dataset2;	/* Dataset ID */
    hsize_t dims1[SPACE1_RANK] =
    {SPACE1_DIM1,SPACE1_DIM2};	/* dataspace dim sizes */
    DATATYPE data_array1[SPACE1_DIM1][SPACE1_DIM2];	/* data buffer */

    hsize_t start[SPACE1_RANK];			/* for hyperslab setting */
    hsize_t count[SPACE1_RANK], stride[SPACE1_RANK];	/* for hyperslab setting */

    herr_t ret;         	/* Generic return value */

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;

    if (verbose)
        printf("Independent write test on file %s\n", filename);

    /* -------------------
     * START AN HDF5 FILE
     * -------------------*/
    /* setup file access template with parallel IO access. */
    acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
    assert(acc_tpl1 != FAIL);
    MESG("H5Pcreate access succeed");
    /* set Parallel access with communicator */
    ret = H5Pset_fapl_mpio(acc_tpl1, comm, info);
    assert(ret != FAIL);
    MESG("H5Pset_fapl_mpio succeed");

    /* create the file collectively */
    fid1 = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl1);
    assert(fid1 != FAIL);
    MESG("H5Fcreate succeed");

    /* Release file-access template */
    ret = H5Pclose(acc_tpl1);
    assert(ret != FAIL);


    /* --------------------------
     * Define the dimensions of the overall datasets
     * and the slabs local to the MPI process.
     * ------------------------- */
    /* setup dimensionality object */
    sid1 = H5Screate_simple(SPACE1_RANK, dims1, NULL);
    assert (sid1 != FAIL);
    MESG("H5Screate_simple succeed");


    /* create a dataset collectively */
    dataset1 = H5Dcreate2(fid1, DATASETNAME1, H5T_NATIVE_INT, sid1,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dataset1 != FAIL);
    MESG("H5Dcreate2 succeed");

    /* create another dataset collectively */
    dataset2 = H5Dcreate2(fid1, DATASETNAME2, H5T_NATIVE_INT, sid1,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dataset2 != FAIL);
    MESG("H5Dcreate2 succeed");



    /* set up dimensions of the slab this process accesses */
    start[0] = mpi_rank*SPACE1_DIM1/mpi_size;
    start[1] = 0;
    count[0] = SPACE1_DIM1/mpi_size;
    count[1] = SPACE1_DIM2;
    stride[0] = 1;
    stride[1] =1;
    if (verbose)
        printf("start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu\n",
                (unsigned long)start[0], (unsigned long)start[1],
                (unsigned long)count[0], (unsigned long)count[1],
                (unsigned long)(count[0]*count[1]));

    /* put some trivial data in the data_array */
    dataset_fill(start, count, stride, &data_array1[0][0]);
    MESG("data_array initialized");

    /* create a file dataspace independently */
    file_dataspace = H5Dget_space (dataset1);
    assert(file_dataspace != FAIL);
    MESG("H5Dget_space succeed");
    ret=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, stride,
            count, NULL);
    assert(ret != FAIL);
    MESG("H5Sset_hyperslab succeed");

    /* create a memory dataspace independently */
    mem_dataspace = H5Screate_simple (SPACE1_RANK, count, NULL);
    assert (mem_dataspace != FAIL);

    /* write data independently */
    ret = H5Dwrite(dataset1, H5T_NATIVE_INT, mem_dataspace, file_dataspace,
            H5P_DEFAULT, data_array1);
    assert(ret != FAIL);
    MESG("H5Dwrite succeed");

    /* write data independently */
    ret = H5Dwrite(dataset2, H5T_NATIVE_INT, mem_dataspace, file_dataspace,
            H5P_DEFAULT, data_array1);
    assert(ret != FAIL);
    MESG("H5Dwrite succeed");

    /* release dataspace ID */
    H5Sclose(file_dataspace);

    /* close dataset collectively */
    ret=H5Dclose(dataset1);
    assert(ret != FAIL);
    MESG("H5Dclose1 succeed");
    ret=H5Dclose(dataset2);
    assert(ret != FAIL);
    MESG("H5Dclose2 succeed");

    /* release all IDs created */
    H5Sclose(sid1);

    /* close the file collectively */
    H5Fclose(fid1);
}

/* Example of using the parallel HDF5 library to read a dataset */
    void
phdf5readInd(char *filename)
{
    hid_t fid1;			/* HDF5 file IDs */
    hid_t acc_tpl1;		/* File access templates */
    hid_t file_dataspace;	/* File dataspace ID */
    hid_t mem_dataspace;	/* memory dataspace ID */
    hid_t dataset1, dataset2;	/* Dataset ID */
    DATATYPE data_array1[SPACE1_DIM1][SPACE1_DIM2];	/* data buffer */
    DATATYPE data_origin1[SPACE1_DIM1][SPACE1_DIM2];	/* expected data buffer */

    hsize_t start[SPACE1_RANK];			/* for hyperslab setting */
    hsize_t count[SPACE1_RANK], stride[SPACE1_RANK];	/* for hyperslab setting */

    herr_t ret;         	/* Generic return value */

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;

    if (verbose)
        printf("Independent read test on file %s\n", filename);

    /* setup file access template */
    acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
    assert(acc_tpl1 != FAIL);
    /* set Parallel access with communicator */
    ret = H5Pset_fapl_mpio(acc_tpl1, comm, info);
    assert(ret != FAIL);


    /* open the file collectively */
    fid1=H5Fopen(filename,H5F_ACC_RDWR,acc_tpl1);
    assert(fid1 != FAIL);

    /* Release file-access template */
    ret=H5Pclose(acc_tpl1);
    assert(ret != FAIL);

    /* open the dataset1 collectively */
    dataset1 = H5Dopen2(fid1, DATASETNAME1, H5P_DEFAULT);
    assert(dataset1 != FAIL);

    /* open another dataset collectively */
    dataset2 = H5Dopen2(fid1, DATASETNAME1, H5P_DEFAULT);
    assert(dataset2 != FAIL);


    /* set up dimensions of the slab this process accesses */
    start[0] = mpi_rank*SPACE1_DIM1/mpi_size;
    start[1] = 0;
    count[0] = SPACE1_DIM1/mpi_size;
    count[1] = SPACE1_DIM2;
    stride[0] = 1;
    stride[1] =1;
    if (verbose)
        printf("start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu\n",
                (unsigned long)start[0], (unsigned long)start[1],
                (unsigned long)count[0], (unsigned long)count[1],
                (unsigned long)(count[0]*count[1]));

    /* create a file dataspace independently */
    file_dataspace = H5Dget_space (dataset1);
    assert(file_dataspace != FAIL);
    ret=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, stride,
            count, NULL);
    assert(ret != FAIL);

    /* create a memory dataspace independently */
    mem_dataspace = H5Screate_simple (SPACE1_RANK, count, NULL);
    assert (mem_dataspace != FAIL);

    /* fill dataset with test data */
    dataset_fill(start, count, stride, &data_origin1[0][0]);

    /* read data independently */
    ret = H5Dread(dataset1, H5T_NATIVE_INT, mem_dataspace, file_dataspace,
            H5P_DEFAULT, data_array1);
    assert(ret != FAIL);

    /* verify the read data with original expected data */
    ret = dataset_vrfy(start, count, stride, &data_array1[0][0], &data_origin1[0][0]);
    assert(ret != FAIL);

    /* read data independently */
    ret = H5Dread(dataset2, H5T_NATIVE_INT, mem_dataspace, file_dataspace,
            H5P_DEFAULT, data_array1);
    assert(ret != FAIL);

    /* verify the read data with original expected data */
    ret = dataset_vrfy(start, count, stride, &data_array1[0][0], &data_origin1[0][0]);
    assert(ret == 0);

    /* close dataset collectively */
    ret=H5Dclose(dataset1);
    assert(ret != FAIL);
    ret=H5Dclose(dataset2);
    assert(ret != FAIL);

    /* release all IDs created */
    H5Sclose(file_dataspace);

    /* close the file collectively */
    H5Fclose(fid1);
}


/*
 * Example of using the parallel HDF5 library to create two datasets
 * in one HDF5 file with collective parallel access support.
 * The Datasets are of sizes (number-of-mpi-processes x DIM1) x DIM2.
 * Each process controls only a slab of size DIM1 x DIM2 within each
 * dataset. [Note: not so yet.  Datasets are of sizes DIM1xDIM2 and
 * each process controls a hyperslab within.]
 */

    void
phdf5writeAll(char *filename,result_pointers_diagnol* result_data)
{
    hid_t fid1;			/* HDF5 file IDs */
    hid_t acc_tpl1;		/* File access templates */
    hid_t xfer_plist;		/* Dataset transfer properties list */
    hid_t sid1;   		/* Dataspace ID */
    hid_t file_dataspace;	/* File dataspace ID */
    hid_t mem_dataspace;	/* memory dataspace ID */
    hid_t dataset_wu, dataset_normal, dataset_sarika, dataset_generalised, dataset_cosine;	/* Dataset ID */
    //hsize_t dims1[SPACE1_RANK] ={SPACE1_DIM1,SPACE1_DIM2};	/* dataspace dim sizes */
    hsize_t dims1[SPACE1_RANK] ={0};	/* dataspace dim sizes */
    DATATYPE data_array1[SPACE1_DIM1][SPACE1_DIM2];	/* data buffer */

    hsize_t start[SPACE1_RANK];			/* for hyperslab setting */
    hsize_t count[SPACE1_RANK], stride[SPACE1_RANK];	/* for hyperslab setting */

    herr_t ret;         	/* Generic return value */

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;

    if (debug_info) {
        if ( debug_mpi_rank == mpi_rank) {
            printf("mpi_rk[%d],result_data->vec_dim=%d\n",mpi_rank,result_data->vec_dim);
            //printf("mpi_rk[%d]:\n",mpi_rank);
            printf("mpi_rk[%d]:%s\n",mpi_rank,result_data->wu_name);
            print_matrix_1d_real(result_data->wu,result_data->vec_dim, "%7.3f");
        }
    }


    if (verbose)
        printf("Collective write test on file %s\n", filename);

    /* -------------------
     * START AN HDF5 FILE
     * -------------------*/
    /* setup file access template with parallel IO access. */
    acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
    assert(acc_tpl1 != FAIL);
    MESG("H5Pcreate access succeed");
    /* set Parallel access with communicator */
    ret = H5Pset_fapl_mpio(acc_tpl1, comm, info);
    assert(ret != FAIL);
    MESG("H5Pset_fapl_mpio succeed");

    /* create the file collectively */
    fid1=H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,acc_tpl1);
    assert(fid1 != FAIL);
    MESG("H5Fcreate succeed");

    /* Release file-access template */
    ret=H5Pclose(acc_tpl1);
    assert(ret != FAIL);


    /* --------------------------
     * Define the dimensions of the overall datasets
     * and create the dataset
     * ------------------------- */
    /* setup dimensionality object */
    dims1[0] = 1; /* this is the total size of the overall dataset */
    dims1[1] = result_data->total_lines*(result_data->total_lines-1)/2;
    sid1 = H5Screate_simple (SPACE1_RANK, dims1, NULL);
    assert (sid1 != FAIL);
    MESG("H5Screate_simple succeed");


    /* create a dataset collectively */
    dataset_wu = H5Dcreate2(fid1, result_data->wu_name, H5T_NATIVE_FLOAT, sid1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dataset_wu != FAIL);
    MESG("H5Dcreate2 succeed");

    /*
     * Set up dimensions of the slab this process accesses.
     */

    /* Dataset1: each process takes a block of rows. */
    //stride = NULL;
    start[0] = 0;
    start[1] = result_data->start_loc;
    if (verbose)
        printf("rk[%d]:start[1]=%d\n",mpi_rank,result_data->start_loc);
    count[0] = 1;
    count[1] = result_data->vec_dim;
    //slab_set(start, count, stride, BYROW);
    if (verbose)
        printf("start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu\n",
                (unsigned long)start[0], (unsigned long)start[1],
                (unsigned long)count[0], (unsigned long)count[1],
                (unsigned long)(count[0]*count[1]));

    /* create a file dataspace independently */
    file_dataspace = H5Dget_space (dataset_wu);
    assert(file_dataspace != FAIL);
    MESG("H5Dget_space succeed");
    //ret=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, stride,
    //    count, NULL);
    ret=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, NULL,
            count, NULL);
    assert(ret != FAIL);
    MESG("H5Sset_hyperslab succeed");

    /* create a memory dataspace independently */
    mem_dataspace = H5Screate_simple (SPACE1_RANK, count, NULL);
    assert (mem_dataspace != FAIL);

    /* fill the local slab with some trivial data */
    //dataset_fill(start, count, stride, &data_array1[0][0]);
    //MESG("data_array initialized");
    if (verbose){
        MESG("data_array created");
        dataset_print(start, count, stride, &data_array1[0][0]);
    }

    /* set up the collective transfer properties list */
    xfer_plist = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist != FAIL);
    ret=H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    assert(ret != FAIL);
    MESG("H5Pcreate xfer succeed");

    /* Create a dataset attribute. */
    /*
     * Create dataspace for the first attribute.
     */
    hsize_t attr_data[]={result_data->total_lines,mpi_size};
    hsize_t adim[] = {sizeof(attr_data)/sizeof(hsize_t)};
    hid_t aid1 = H5Screate(H5S_SIMPLE);
    ret  = H5Sset_extent_simple(aid1, 1, adim, NULL);
    assert(ret != FAIL);
    /*
     * Create array attribute.
     */
    //attr1 = H5Acreate2(dataset, ANAME, H5T_NATIVE_FLOAT, aid1, H5P_DEFAULT, H5P_DEFAULT);
    hid_t root_grp = H5Gopen1(fid1,"/");
    //hid_t attribute_id = H5Acreate2 (root_grp, "MatrixInfo", H5T_STD_I32BE, aid1, 
    //        H5P_DEFAULT, H5P_DEFAULT);
    hid_t attribute_id = H5Acreate2 (root_grp, "MatrixInfo", H5T_STD_U64LE, aid1, 
            H5P_DEFAULT, H5P_DEFAULT);
    /* Write the attribute data. */
    ret = H5Awrite(attribute_id, H5T_NATIVE_ULLONG, attr_data);
    assert(ret != FAIL);
    MESG("H5Awrite succeed");
    /* Close the attribute. */
    ret = H5Aclose(attribute_id);
    assert(ret != FAIL);
    H5Gclose(root_grp);


    /* write data collectively */
    //ret = H5Dwrite(dataset_wu, H5T_NATIVE_INT, mem_dataspace, file_dataspace,
    //    xfer_plist, data_array1);
    ret = H5Dwrite(dataset_wu, H5T_NATIVE_FLOAT, mem_dataspace, file_dataspace,
            xfer_plist, result_data->wu);
    assert(ret != FAIL);
    MESG("H5Dwrite succeed");

    dataset_normal = H5Dcreate2(fid1, result_data->normal_name, H5T_NATIVE_FLOAT, sid1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dataset_normal != FAIL);
    MESG("H5Dcreate2 succeed");
    ret = H5Dwrite(dataset_normal, H5T_NATIVE_FLOAT, mem_dataspace, file_dataspace,
            xfer_plist, result_data->normal);
    assert(ret != FAIL);
    MESG("H5Dwrite succeed");



    dataset_sarika = H5Dcreate2(fid1, result_data->sarika_name, H5T_NATIVE_FLOAT, sid1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dataset_sarika != FAIL);
    MESG("H5Dcreate2 succeed");
    ret = H5Dwrite(dataset_sarika, H5T_NATIVE_FLOAT, mem_dataspace, file_dataspace,
            xfer_plist, result_data->sarika);
    assert(ret != FAIL);
    MESG("H5Dwrite succeed");

    dataset_generalised = H5Dcreate2(fid1, result_data->generalised_name, H5T_NATIVE_FLOAT, sid1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dataset_generalised != FAIL);
    MESG("H5Dcreate2 succeed");
    ret = H5Dwrite(dataset_generalised, H5T_NATIVE_FLOAT, mem_dataspace, file_dataspace,
            xfer_plist, result_data->generalised);
    assert(ret != FAIL);
    MESG("H5Dwrite succeed");

    dataset_cosine = H5Dcreate2(fid1, result_data->cosine_name, H5T_NATIVE_FLOAT, sid1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dataset_cosine != FAIL);
    MESG("H5Dcreate2 succeed");
    ret = H5Dwrite(dataset_cosine, H5T_NATIVE_FLOAT, mem_dataspace, file_dataspace,
            xfer_plist, result_data->cosine);
    assert(ret != FAIL);
    MESG("H5Dwrite succeed");


    /* release all temporary handles. */
    /* Could have used them for dataset2 but it is cleaner */
    /* to create them again.*/
    H5Sclose(file_dataspace);
    H5Sclose(mem_dataspace);
    H5Pclose(xfer_plist);

    //    /*
    //     * All writes completed.  Close datasets collectively
    //     */
    ret=H5Dclose(dataset_wu);
    assert(ret != FAIL);
    ret=H5Dclose(dataset_normal);
    assert(ret != FAIL);
    ret=H5Dclose(dataset_sarika);
    assert(ret != FAIL);
    ret=H5Dclose(dataset_generalised);
    assert(ret != FAIL);
    ret=H5Dclose(dataset_cosine);
    assert(ret != FAIL);

    /* write the start_loc info
    /* setup dimensionality object */
    dims1[0] = 6; /* this is the total size of the overall dataset */
    dims1[1] = mpi_size;
    sid1 = H5Screate_simple (SPACE1_RANK, dims1, NULL);
    assert (sid1 != FAIL);
    MESG("H5Screate_simple succeed");

    /* create a dataset collectively */
    //hid_t dataset_sl_cnt = H5Dcreate2(fid1, "start_loc", H5T_NATIVE_INT, sid1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t dataset_sl_cnt = H5Dcreate2(fid1, "start_loc", H5T_STD_U64LE, sid1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dataset_sl_cnt != FAIL);
    MESG("H5Dcreate2 succeed");
    start[0] = 0;
    start[1] = mpi_rank;
    count[0] = 6;
    count[1] = 1;

    /* create a file dataspace independently */
    file_dataspace = H5Dget_space (dataset_sl_cnt);
    assert(file_dataspace != FAIL);
    MESG("H5Dget_space succeed");
    //ret=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, stride,
    //    count, NULL);
    ret=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, NULL,
            count, NULL);
    assert(ret != FAIL);
    MESG("H5Sset_hyperslab succeed");

    /* create a memory dataspace independently */
    mem_dataspace = H5Screate_simple (SPACE1_RANK, count, NULL);
    assert (mem_dataspace != FAIL);
    /* set up the collective transfer properties list */
    xfer_plist = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist != FAIL);
    ret=H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    assert(ret != FAIL);
    MESG("H5Pcreate xfer succeed");

    //dataset_normal = H5Dcreate2(fid1, result_data->normal_name, H5T_NATIVE_FLOAT, sid1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //assert(dataset_normal != FAIL);
    //MESG("H5Dcreate2 succeed");
    integer sl_cnt[6][1] = {{result_data->start_loc},{result_data->vec_dim},
        {result_data->chunk_start_a},
        {result_data->chunk_start_b},
        {result_data->chunk_count_a},
        {result_data->chunk_count_b}};
    //ret = H5Dwrite(dataset_sl_cnt, H5T_NATIVE_INT, mem_dataspace, file_dataspace,
    //        xfer_plist, sl_cnt);
    ret = H5Dwrite(dataset_sl_cnt, H5T_NATIVE_ULLONG, mem_dataspace, file_dataspace,
            xfer_plist, sl_cnt);
    assert(ret != FAIL);
    MESG("H5Dwrite succeed");

    ret=H5Dclose(dataset_sl_cnt);
    assert(ret != FAIL);

    H5Sclose(file_dataspace);
    H5Sclose(mem_dataspace);
    H5Pclose(xfer_plist);

    /* release all IDs created */
    H5Sclose(sid1);

    /* close the file collectively */
    H5Fclose(fid1);
}

/*
 * Example of using the parallel HDF5 library to read two datasets
 * in one HDF5 file with collective parallel access support.
 * The Datasets are of sizes (number-of-mpi-processes x DIM1) x DIM2.
 * Each process controls only a slab of size DIM1 x DIM2 within each
 * dataset. [Note: not so yet.  Datasets are of sizes DIM1xDIM2 and
 * each process controls a hyperslab within.]
 */

    void
phdf5readAll(char *filename)
{
    hid_t fid1;			/* HDF5 file IDs */
    hid_t acc_tpl1;		/* File access templates */
    hid_t xfer_plist;		/* Dataset transfer properties list */
    hid_t file_dataspace;	/* File dataspace ID */
    hid_t mem_dataspace;	/* memory dataspace ID */

    hsize_t     dims_out[2];

    hid_t dataset1, dataset2;	/* Dataset ID */
    //DATATYPE data_array_a[SPACE1_DIM1][SPACE1_DIM2];	/* data buffer */
    //unsigned long space_dim_a0=0, space_dim_a1=0;
    //unsigned long space_dim_b0=0, space_dim_b1=0;
    hsize_t space_dim_a0=0, space_dim_a1=0;
    hsize_t space_dim_b0=0, space_dim_b1=0;
    //DATATYPE data_array_a[space_dim_a1][space_dim2];	/* data buffer */
    //DATATYPE **data_array_a=NULL;	/* data buffer */
    DATATYPE data_origin1[SPACE1_DIM1][SPACE1_DIM2];	/* expected data buffer */

    hsize_t start[SPACE1_RANK];			/* for hyperslab setting */
    hsize_t count[SPACE1_RANK], stride[SPACE1_RANK];	/* for hyperslab setting */

    hsize_t start_part_a[SPACE1_RANK], start_part_b[SPACE1_RANK];			/* for hyperslab setting */
    hsize_t count_part_a[SPACE1_RANK], count_part_b[SPACE1_RANK];

    herr_t ret;         	/* Generic return value */
    integer i,j,status_n, fs_rank;
    integer average_lines,remainder_lines, num_data_chunks;
    double t1, t2; 


    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;

    MPI_Barrier(comm);
    t1 = MPI_Wtime(); 

    if (verbose)
        printf("Collective read test on file %s\n", filename);

    /* -------------------
     * OPEN AN HDF5 FILE
     * -------------------*/
    /* setup file access template with parallel IO access. */
    acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
    assert(acc_tpl1 != FAIL);
    MESG("H5Pcreate access succeed");
    /* set Parallel access with communicator */
    ret = H5Pset_fapl_mpio(acc_tpl1, comm, info);
    assert(ret != FAIL);
    MESG("H5Pset_fapl_mpio succeed");

    /* open the file collectively */
    fid1=H5Fopen(filename,H5F_ACC_RDWR,acc_tpl1);
    assert(fid1 != FAIL);
    MESG("H5Fopen succeed");

    /* Release file-access template */
    ret=H5Pclose(acc_tpl1);
    assert(ret != FAIL);


    /* --------------------------
     * Open the datasets in it
     * ------------------------- */
    /* open the dataset1 collectively */
    dataset1 = H5Dopen2(fid1, DATASETNAME1, H5P_DEFAULT);
    assert(dataset1 != FAIL);
    MESG("H5Dopen2 succeed");

    /* open another dataset collectively */
    //   dataset2 = H5Dopen2(fid1, DATASETNAME1, H5P_DEFAULT);
    //   assert(dataset2 != FAIL);
    //   MESG("H5Dopen2 2 succeed");

    /*
     * Set up dimensions of the slab this process accesses.
     */

    /* Dataset1: each process takes a block of rows. */
    //slab_set(start, count, stride, BYROW);

    /* create a file dataspace independently */
    file_dataspace = H5Dget_space (dataset1);
    fs_rank      = H5Sget_simple_extent_ndims (file_dataspace);
    status_n  = H5Sget_simple_extent_dims (file_dataspace, dims_out, NULL);
    if (verbose)
        printf("mat_rk=%d,dims_out[0]=%lu,dims_out[1]=%lu\n", 
                fs_rank,
                (unsigned long)dims_out[0],(unsigned long)dims_out[1]);
    assert(file_dataspace != FAIL);
    MESG("H5Dget_space succeed");

    // Calculate number of data chunks
    num_data_chunks = (int)sqrt(mpi_size*2);
    if (num_data_chunks*num_data_chunks != mpi_size*2){
        if (0 == mpi_rank) {
            printf("Error: np*2 is not a square number! Exiting...\n");
        }
        MPI_Finalize();
        exit (0);
    }
    if (verbose)
        printf("mpi_size=%d,ndchunks=%d\n",mpi_size,num_data_chunks);
    // assign chunk combinations to each mpi_rank
    //int n=5,r=2; 
    integer **cmbs=NULL;
    integer num_cmbs = 0;
    cmbs=combination_util(num_data_chunks,&num_cmbs); 
    //printf("num_cmb=%d,num_cmbs1=%d\n",num_cmbs,num_cmbs1);
    //print_matrix(cmbs, num_cmbs, 2, "%3d");
    integer mpi_rk_chunk0, mpi_rk_chunk1;
    if (mpi_rank < num_cmbs) { // an off-diagnal full block
        //printf("mpi_rank=%d\n",mpi_rank);
        mpi_rk_chunk0 = cmbs[mpi_rank][0];
        mpi_rk_chunk1 = cmbs[mpi_rank][1];
    }
    else { // two diagnal triangles
        mpi_rk_chunk0 = (mpi_rank-num_cmbs)*2;
        mpi_rk_chunk1 = mpi_rk_chunk0 + 1;
    }
    if (verbose)
        printf("mrk[%d]:mpi_rk_chunk0=%d, mpi_rk_chunk1=%d\n",mpi_rank,mpi_rk_chunk0, mpi_rk_chunk1);
    free_dynamic_2d_array_integer(cmbs);

    /* now calculate start[0] for each rank*/
    average_lines = dims_out[0] / num_data_chunks;
    remainder_lines = dims_out[0] % num_data_chunks;
    if (verbose)
        printf("average_lines=%d, remainder_lines=%d\n", average_lines,remainder_lines);

    integer **chunk_start=NULL, **chunk_count=NULL;
    chunk_start = allocate_dynamic_2d_array_integer(num_data_chunks,2);
    chunk_count = allocate_dynamic_2d_array_integer(num_data_chunks,2);
    for (i=0;i<remainder_lines;i++) {
        chunk_start[i][0] = i*(average_lines+1);
        chunk_start[i][1] = 0;
        chunk_count[i][0] = average_lines+1;
        chunk_count[i][1] = dims_out[1];
    }

    for (i=remainder_lines;i<num_data_chunks;i++){
        chunk_start[i][0] = remainder_lines*(average_lines+1)+(i-remainder_lines)*average_lines;
        chunk_start[i][1] = 0;
        chunk_count[i][0] = average_lines;
        chunk_count[i][1] = dims_out[1];
    }

    if (verbose)
        if ( 0==mpi_rank ) {
            print_matrix_integer(chunk_start,num_data_chunks,2,"%3d ");
            print_matrix_integer(chunk_count,num_data_chunks,2,"%3d ");
        }

    // now starting assign two chunks, part_a and part_b to each mpi rank
    // you might directly use chunk_start[mpi_rk_chunk0] and chunk_start[mpi_rk_chunk1]
    // but I believe this more understandable

    // part_a
    start_part_a[0]=chunk_start[mpi_rk_chunk0][0];
    start_part_a[1]=chunk_start[mpi_rk_chunk0][1];
    count_part_a[0]=chunk_count[mpi_rk_chunk0][0];
    count_part_a[1]=chunk_count[mpi_rk_chunk0][1];

    // part b
    start_part_b[0]=chunk_start[mpi_rk_chunk1][0];
    start_part_b[1]=chunk_start[mpi_rk_chunk1][1];
    count_part_b[0]=chunk_count[mpi_rk_chunk1][0];
    count_part_b[1]=chunk_count[mpi_rk_chunk1][1];

    if (verbose)
        printf("start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu\n",
                (unsigned long)start[0], (unsigned long)start[1],
                (unsigned long)count[0], (unsigned long)count[1],
                (unsigned long)(count[0]*count[1]));
    ret=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start_part_a, NULL,
            count_part_a, NULL);
    assert(ret != FAIL);
    MESG("H5Sset_hyperslab succeed");
    //
    ///* create a memory dataspace independently */
    mem_dataspace = H5Screate_simple (SPACE1_RANK, count_part_a, NULL);
    assert (mem_dataspace != FAIL);
    //
    //DATATYPE **data_array_a=NULL;	/* data buffer */
    //DATATYPE **data_array_b=NULL;
    sint **data_array_a=NULL;	/* data buffer */
    sint **data_array_b=NULL;
    //space_dim_a0 = (unsigned long)count_part_a[0];
    //space_dim_a1 = (unsigned long)count_part_a[1];
    //space_dim_b0 = (unsigned long)count_part_b[0];
    //space_dim_b1 = (unsigned long)count_part_b[1];
    space_dim_a0 = (hsize_t)count_part_a[0];
    space_dim_a1 = (hsize_t)count_part_a[1];
    space_dim_b0 = (hsize_t)count_part_b[0];
    space_dim_b1 = (hsize_t)count_part_b[1];
    if (verbose)
    {
        printf("space_dim_a0=%lu,space_dim_a1=%lu\n",space_dim_a0,space_dim_a1);
        printf("space_dim_b0=%lu,space_dim_b1=%lu\n",space_dim_b0,space_dim_b1);
    }

    //data_array_a = allocate_dynamic_2d_array_integer(space_dim_a0,space_dim_a1);
    //data_array_b = allocate_dynamic_2d_array_integer(space_dim_a0,space_dim_a1);
    data_array_a = allocate_dynamic_2d_array_sint(space_dim_a0,space_dim_a1);
    data_array_b = allocate_dynamic_2d_array_sint(space_dim_b0,space_dim_b1);

    ///* set up the collective transfer properties list */
    xfer_plist = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist != FAIL);
    ret=H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    assert(ret != FAIL);
    MESG("H5Pcreate xfer succeed");
    //
    ///* read data collectively */
    ret = H5Dread(dataset1, H5T_NATIVE_SHORT, mem_dataspace, file_dataspace,
            xfer_plist, &data_array_a[0][0]);
    assert(ret != FAIL);
    MESG("H5Dread data_array_a succeed");

    //printf("data_array_a=%p\n",data_array_a);
    //printf("%5d",data_array_a[0][0]);



    ret=H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start_part_b, NULL,
            count_part_b, NULL);
    assert(ret != FAIL);
    MESG("H5Sset_hyperslab data_array_b succeed");

    mem_dataspace = H5Screate_simple (SPACE1_RANK, count_part_b, NULL);
    assert (mem_dataspace != FAIL);

    ret = H5Dread(dataset1, H5T_NATIVE_SHORT, mem_dataspace, file_dataspace,
            xfer_plist, &data_array_b[0][0]);
    assert(ret != FAIL);
    MESG("H5Dread data_array_b succeed");

    // verify results
    //int debug_mpi_rank = 7;
    if (debug_info)
        if (debug_mpi_rank == mpi_rank) 
        {
            printf("debug_mpi_rank=%d\n",debug_mpi_rank);
            for (i=0;i<space_dim_a0;i++)
            {
                printf("mpi_rank[%d]a:",mpi_rank);
                for (j=0;j<space_dim_a1;j++)
                    printf("%4d,",data_array_a[i][j]);
                printf("\n");
            }
            for (i=0;i<space_dim_b0;i++)
            {
                printf("mpi_rank[%d]b:",mpi_rank);
                for (j=0;j<space_dim_b1;j++)
                    printf("%4d,",data_array_b[i][j]);
                printf("\n");
            }
        }

    /* release all temporary handles. */
    /* Could have used them for dataset2 but it is cleaner */
    /* to create them again.*/
    H5Sclose(file_dataspace);
    H5Sclose(mem_dataspace);
    H5Pclose(xfer_plist);

    /*
     * All reads completed.  Close datasets collectively
     */
    ret=H5Dclose(dataset1);
    assert(ret != FAIL);
    MESG("H5Dclose1 succeed");

    /* close the file collectively */
    H5Fclose(fid1);

    MPI_Barrier(comm);
    t2 = MPI_Wtime(); 
    if (0==mpi_rank)
        printf( "Reading time is %.3f\n", t2 - t1 ); 

    MPI_Barrier(comm);
    t1 = MPI_Wtime(); 

    /* now starting on the processing */
    real *part_ab_normal = NULL; 
    real *part_ab_generalised = NULL; 
    real *part_ab_wu = NULL; 
    real *part_ab_sarika = NULL; 
    real *part_ab_cosine = NULL; 
    integer num_cmbs_ab;

    if (mpi_rank < num_cmbs){
        //real **block_cmbs;
        /* process off-diagnol blocks*/

        num_cmbs_ab = space_dim_a0*space_dim_b0;
        part_ab_normal = (real*)malloc(num_cmbs_ab*sizeof(real));
        part_ab_generalised = (real*)malloc(num_cmbs_ab*sizeof(real));
        part_ab_wu = (real*)malloc(num_cmbs_ab*sizeof(real));
        part_ab_sarika = (real*)malloc(num_cmbs_ab*sizeof(real));
        part_ab_cosine = (real*)malloc(num_cmbs_ab*sizeof(real));

        result_pointers_diagnol rp;
        rp.normal = part_ab_normal;
        rp.generalised = part_ab_generalised;
        rp.wu = part_ab_wu;
        rp.sarika = part_ab_sarika;
        rp.cosine = part_ab_cosine;

        calc_coeffs_off_diagnol_block(data_array_a, space_dim_a0, space_dim_a1, 
                data_array_b, space_dim_b0, space_dim_b1,
                &rp);
        if (debug_info)
            if (mpi_rank == debug_mpi_rank) {
                // we can treat as 1D array as well
                //print_matrix_1d_real(rp.normal[0], space_dim_a0*space_dim_b0, "%7.3f ");
                printf("mpi_rank[%d]:normal",mpi_rank);
                print_matrix_1d_real(rp.normal, num_cmbs_ab, "%7.3f ");
                printf("mpi_rank[%d]:generalised",mpi_rank);
                print_matrix_1d_real(rp.generalised, num_cmbs_ab, "%7.3f ");
                printf("mpi_rank[%d]:wu",mpi_rank);
                print_matrix_1d_real(rp.wu, num_cmbs_ab, "%7.3f ");
                printf("mpi_rank[%d]:sarika",mpi_rank);
                print_matrix_1d_real(rp.sarika, num_cmbs_ab, "%7.3f ");
                printf("mpi_rank[%d]:cosine",mpi_rank);
                print_matrix_1d_real(rp.cosine, num_cmbs_ab, "%7.3f ");
            }

        /* writes to hdf5 collectively*/


        //free_dynamic_2d_array_real(rp.normal);
        //free_dynamic_2d_array_real(rp.generalised);
        //free_dynamic_2d_array_real(rp.wu);
        //free_dynamic_2d_array_real(rp.sarika);
        //free_dynamic_2d_array_real(rp.cosine);

    }
    else {
        /* process two diagnol triangles*/

        result_pointers_diagnol rpd_part_a;
        integer num_cmbs_a = space_dim_a0*(space_dim_a0-1)/2;
        integer num_cmbs_b = space_dim_b0*(space_dim_b0-1)/2;
        num_cmbs_ab = num_cmbs_a + num_cmbs_b;

        part_ab_normal = (real*)malloc(num_cmbs_ab*sizeof(real)); 
        part_ab_generalised = (real*)malloc(num_cmbs_ab*sizeof(real)); 
        part_ab_wu = (real*)malloc(num_cmbs_ab*sizeof(real));
        part_ab_sarika = (real*)malloc(num_cmbs_ab*sizeof(real)); 
        part_ab_cosine = (real*)malloc(num_cmbs_ab*sizeof(real)); 

        rpd_part_a.normal = part_ab_normal;
        rpd_part_a.generalised = part_ab_generalised;
        rpd_part_a.wu = part_ab_wu;
        rpd_part_a.sarika = part_ab_sarika;
        rpd_part_a.cosine = part_ab_cosine;

        //calc_coeffs_diagnol_triangle(data_array_a, space_dim_a0, space_dim_a1,&rpd_part_a);
        calc_coeffs_off_diagnol_block(data_array_a, space_dim_a0, space_dim_a1, 
                data_array_a, space_dim_a0, space_dim_a1,
                &rpd_part_a);

        if (debug_info)
            if (mpi_rank == debug_mpi_rank) {
                printf("mpi_rank[%d]:normal_part_a",mpi_rank);
                print_matrix_1d_real(rpd_part_a.normal, num_cmbs_a, "%7.3f ");
                printf("mpi_rank[%d]:generalised_part_a",mpi_rank);
                print_matrix_1d_real(rpd_part_a.generalised, num_cmbs_a, "%7.3f ");
                printf("mpi_rank[%d]:wu_part_a",mpi_rank);
                print_matrix_1d_real(rpd_part_a.wu, num_cmbs_a, "%7.3f ");
                printf("mpi_rank[%d]:sarika_part_a",mpi_rank);
                print_matrix_1d_real(rpd_part_a.sarika, num_cmbs_a, "%7.3f ");
                printf("mpi_rank[%d]:cosine_part_a",mpi_rank);
                print_matrix_1d_real(rpd_part_a.cosine, num_cmbs_a, "%7.3f ");
            }

        result_pointers_diagnol rpd_part_b;
        // second section of the triangle results
        rpd_part_b.normal = part_ab_normal + num_cmbs_a;
        rpd_part_b.generalised = part_ab_generalised + num_cmbs_a;
        rpd_part_b.wu = part_ab_wu + num_cmbs_a;
        rpd_part_b.sarika = part_ab_sarika + num_cmbs_a;
        rpd_part_b.cosine = part_ab_cosine + num_cmbs_a;

        //calc_coeffs_diagnol_triangle(data_array_b, space_dim_b0, space_dim_b1,&rpd_part_b);
        calc_coeffs_off_diagnol_block(data_array_b, space_dim_b0, space_dim_b1, 
                data_array_b, space_dim_b0, space_dim_b1,
                &rpd_part_b);

        if (debug_info)
            if (mpi_rank == debug_mpi_rank) {
                printf("mpi_rank[%d]:normal_part_b",mpi_rank);
                print_matrix_1d_real(rpd_part_b.normal, num_cmbs_b, "%7.3f ");
                printf("mpi_rank[%d]:generalised_part_b",mpi_rank);
                print_matrix_1d_real(rpd_part_b.generalised, num_cmbs_b, "%7.3f ");
                printf("mpi_rank[%d]:wu_part_b",mpi_rank);
                print_matrix_1d_real(rpd_part_b.wu, num_cmbs_b, "%7.3f ");
                printf("mpi_rank[%d]:sarika_part_b",mpi_rank);
                print_matrix_1d_real(rpd_part_b.sarika, num_cmbs_b, "%7.3f ");
                printf("mpi_rank[%d]:cosine_part_b",mpi_rank);
                print_matrix_1d_real(rpd_part_b.cosine, num_cmbs_b, "%7.3f ");
            }
    }

    /* now starting to work on the writing */
    //write_hdf5(result_pointers_diagnol,num_cmbs_ab,mpi_size,mpi_rank);
    result_pointers_diagnol rpd_h5;
    rpd_h5.normal = part_ab_normal;
    rpd_h5.normal_name = "normal";
    rpd_h5.generalised = part_ab_generalised;
    rpd_h5.generalised_name = "generalised";
    rpd_h5.wu = part_ab_wu;
    rpd_h5.wu_name = "wu";
    rpd_h5.sarika = part_ab_sarika;
    rpd_h5.sarika_name = "sarika";
    rpd_h5.cosine = part_ab_cosine;
    rpd_h5.cosine_name = "cosine";
    rpd_h5.vec_dim = num_cmbs_ab;

    rpd_h5.total_lines = dims_out[0]; // total number of lines

    // get the start location of each rank for write
    integer *rk_vec_dim_all=(integer *)malloc(sizeof(integer)*mpi_size);
    integer *rk_start=(integer *)malloc(sizeof(integer)*mpi_size);
    //MPI_Allgather(&num_cmbs_ab,1,MPI_INT,rk_vec_dim_all,1,MPI_INT,comm);
    MPI_Allgather(&num_cmbs_ab,1,MPI_UNSIGNED_LONG_LONG,rk_vec_dim_all,1,MPI_UNSIGNED_LONG_LONG,comm);
    integer cur_loc = 0;
    for (i=0;i<mpi_size;i++) {
        rk_start[i] = cur_loc;
        cur_loc += rk_vec_dim_all[i];
    }
    rpd_h5.start_loc = rk_start[mpi_rank];
    rpd_h5.chunk_start_a = start_part_a[0];
    rpd_h5.chunk_start_b = start_part_b[0];
    rpd_h5.chunk_count_a = count_part_a[0];
    rpd_h5.chunk_count_b = count_part_b[0];

    MPI_Barrier(comm);
    t2 = MPI_Wtime(); 
    if (0==mpi_rank)
        printf( "Elapsed time is %.3f\n", t2 - t1 ); 


    MPI_Barrier(comm);
    t1 = MPI_Wtime(); 
    //phdf5writeAll("res_all.h5",&rpd_h5);
    phdf5writeAll(testfiles[1],&rpd_h5);

    MPI_Barrier(comm);
    t2 = MPI_Wtime(); 
    if (0==mpi_rank)
        printf( "Writing time is %.3f\n", t2 - t1 ); 

    free(rk_vec_dim_all);
    free(rk_start);

    free(part_ab_normal);
    free(part_ab_generalised);
    free(part_ab_wu);
    free(part_ab_sarika);
    free(part_ab_cosine);

    //free_dynamic_2d_array_integer(data_array_a);
    //free_dynamic_2d_array_integer(data_array_b);
    free_dynamic_2d_array_sint(data_array_a);
    free_dynamic_2d_array_sint(data_array_b);
    free_dynamic_2d_array_integer(chunk_start);
    free_dynamic_2d_array_integer(chunk_count);

}

/*
 * test file access by communicator besides COMM_WORLD.
 * Split COMM_WORLD into two, one (even_comm) contains the original
 * processes of even ranks.  The other (odd_comm) contains the original
 * processes of odd ranks.  Processes in even_comm creates a file, then
 * cloose it, using even_comm.  Processes in old_comm just do a barrier
 * using odd_comm.  Then they all do a barrier using COMM_WORLD.
 * If the file creation and cloose does not do correct collective action
 * according to the communicator argument, the processes will freeze up
 * sooner or later due to barrier mixed up.
 */
    void
test_split_comm_access(char filenames[][PATH_MAX])
{
    MPI_Comm comm;
    MPI_Info info = MPI_INFO_NULL;
    int color, mrc;
    int newrank, newprocs;
    hid_t fid;			/* file IDs */
    hid_t acc_tpl;		/* File access properties */
    herr_t ret;			/* generic return value */

    if (verbose)
        printf("Independent write test on file %s %s\n",
                filenames[0], filenames[1]);

    color = mpi_rank%2;
    mrc = MPI_Comm_split (MPI_COMM_WORLD, color, mpi_rank, &comm);
    assert(mrc==MPI_SUCCESS);
    MPI_Comm_size(comm,&newprocs);
    MPI_Comm_rank(comm,&newrank);

    if (color){
        /* odd-rank processes */
        mrc = MPI_Barrier(comm);
        assert(mrc==MPI_SUCCESS);
    }else{
        /* even-rank processes */
        /* setup file access template */
        acc_tpl = H5Pcreate (H5P_FILE_ACCESS);
        assert(acc_tpl != FAIL);

        /* set Parallel access with communicator */
        ret = H5Pset_fapl_mpio(acc_tpl, comm, info);
        assert(ret != FAIL);

        /* create the file collectively */
        fid=H5Fcreate(filenames[color],H5F_ACC_TRUNC,H5P_DEFAULT,acc_tpl);
        assert(fid != FAIL);
        MESG("H5Fcreate succeed");

        /* Release file-access template */
        ret=H5Pclose(acc_tpl);
        assert(ret != FAIL);

        ret=H5Fclose(fid);
        assert(ret != FAIL);
    }
    if (mpi_rank == 0){
        mrc = MPI_File_delete(filenames[color], info);
        assert(mrc==MPI_SUCCESS);
    }
}

/*
 * Show command usage
 */
    void
usage(void)
{
    printf("Usage: testphdf5 [-f <prefix>] [-r] [-w] [-v]\n");
    printf("\t-f\tfile prefix for parallel test files.\n");
    printf("\t  \te.g. pfs:/PFS/myname\n");
    printf("\t  \tcan be set via $" PARAPREFIX ".\n");
    printf("\t  \tDefault is current directory.\n");
    printf("\t-c\tno cleanup\n");
    printf("\t-r\tno read\n");
    printf("\t-w\tno write\n");
    printf("\t-v\tverbose on\n");
    printf("\tdefault do write then read\n");
    printf("\n");
}


/*
 * compose the test filename with the prefix supplied.
 * return code: 0 if no error
 *              1 otherwise.
 */
    int
mkfilenames(char *prefix)
{
    int i, n;
    size_t strsize;

    /* filename will be prefix/ParaEgN.h5 where N is 0 to 9. */
    /* So, string must be big enough to hold the prefix, / and 10 more chars */
    /* and the terminating null. */
    strsize = strlen(prefix) + 12;
    if (strsize > PATH_MAX){
        printf("File prefix too long;  Use a short path name.\n");
        return(1);
    }
    n = sizeof(testfiles)/sizeof(testfiles[0]);
    if (n > 9){
        printf("Warning: Too many entries in testfiles. "
                "Need to adjust the code to accommodate the large size.\n");
    }
    //for (i=0; i<n; i++){
    //    /*sprintf(testfiles[i], "%s/ParaEg%d.h5", prefix, i);*/
    //    sprintf(testfiles[i], "%s", prefix);
    //}
    sprintf(testfiles[0], "%s", prefix);
    sprintf(testfiles[1], "%s.res_all.h5", prefix);
    return(0);

}


/*
 * parse the command line options
 */
int
parse_options(int argc, char **argv){
    int i, n;

    /* initialize testfiles to nulls */
    n = sizeof(testfiles)/sizeof(testfiles[0]);
    for (i=0; i<n; i++){
        testfiles[i][0] = '\0';
    }

    while (--argc){
        if (**(++argv) != '-'){
            break;
        }else{
            switch(*(*argv+1)){
                case 'f':   ++argv;
                            if (--argc < 1){
                                usage();
                                nerrors++;
                                return(1);
                            }
                            if (mkfilenames(*argv)){
                                nerrors++;
                                return(1);
                            }
                            break;
                case 'c':   docleanup = 0;	/* no cleanup */
                            break;
                case 'r':   doread = 0;
                            break;
                case 'w':   dowrite = 0;
                            break;
                case 'v':   verbose = 1;
                            break;
                case 'g':   debug_info = 1;
                            break;
                case 'd':   
                            ++argv,--argc;
                            debug_mpi_rank = atoi(*argv);

                            break;
                default:    usage();
                            nerrors++;
                            return(1);
            }
        }
    }

    /* check the file prefix */
    if (testfiles[0][0] == '\0'){
        /* try get it from environment variable HDF5_PARAPREFIX */
        char *env;
        char *env_default = ".";	/* default to current directory */
        if ((env=getenv(PARAPREFIX))==NULL){
            env = env_default;
        }
        mkfilenames(env);
    }
    return(0);
}


/*
 * cleanup test files created
 */
    void
cleanup(void)
{
    int i, n;

    n = sizeof(testfiles)/sizeof(testfiles[0]);
    for (i=0; i<n; i++){
        MPI_File_delete(testfiles[i], MPI_INFO_NULL);
    }
}


/* Main Program */
    int
main(int argc, char **argv)
{
    int mpi_namelen;
    char mpi_name[MPI_MAX_PROCESSOR_NAME];
    int i, n;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
    MPI_Get_processor_name(mpi_name,&mpi_namelen);
    /* Make sure datasets can be divided into equal chunks by the processes */
    //if ((SPACE1_DIM1 % mpi_size) || (SPACE1_DIM2 % mpi_size)){
    //printf("DIM1(%d) and DIM2(%d) must be multiples of processes (%d)\n",
    //    SPACE1_DIM1, SPACE1_DIM2, mpi_size);
    //nerrors++;
    //goto finish;
    //}

//#ifdef _OPENACC
//    acc_init(acc_device_nvidia);                                 // OpenACC call
//    ngpus = acc_get_num_devices(acc_device_nvidia);  // #GPUs
//    int dev_id = myrank % num_dev;         
//    acc_set_device_num(dev_id,acc_device_nvidia); // assign GPU to one MPI process
//
//    cout << "MPI process " << myrank << "  is assigned to GPU " << dev_id << "\n";
//#endif

#ifdef _OPENACC
    acc_init(acc_device_nvidia);                                 // OpenACC call
    ngpus = acc_get_num_devices( acc_device_nvidia ); 
    if( ngpus ){
        gpunum = mpi_rank % ngpus;
        acc_set_device_num( gpunum, acc_device_nvidia ); 
    }
    else {
        acc_set_device_type( acc_device_host ); 
    }
    printf("MPI process %d is on GPU %d\n", mpi_rank, gpunum);
#endif

    if (parse_options(argc, argv) != 0)
        goto finish;

    /* show test file names */
    if (mpi_rank == 0){
        n = sizeof(testfiles)/sizeof(testfiles[0]);
        //printf("Parallel test files are:\n");
        //for (i=0; i<n; i++){
        //    printf("   %s\n", testfiles[i]);
        //}
        printf("Input file is: %s\n", testfiles[0]);
        printf("Output file is: %s\n", testfiles[1]);
    }

    if (dowrite){
        MPI_BANNER("testing PHDF5 dataset using split communicators...");
        test_split_comm_access(testfiles);
        MPI_BANNER("testing PHDF5 dataset independent write...");
        phdf5writeInd(testfiles[0]);
        MPI_BANNER("testing PHDF5 dataset collective write...");
        phdf5writeAll(testfiles[1],NULL);
    }
    if (doread){
        //MPI_BANNER("testing PHDF5 dataset independent read...");
        //phdf5readInd(testfiles[0]);
        MPI_BANNER("testing PHDF5 dataset collective read...");
        phdf5readAll(testfiles[0]);
        //phdf5readAll(testfiles[1]);
        //phdf5readAll("ta.h5");
        //phdf5readAll("sample_hsp70_actin.h5");
        //phdf5readAll("sample_a-b_mix_2.h5");
        //phdf5readAll("sample_protease_mix_1.h5");
    }

    if (!(dowrite || doread)){
        usage();
        nerrors++;
    }

finish:
    if (mpi_rank == 0){		/* only process 0 reports */
        if (nerrors)
            printf("***PHDF5 tests detected %d errors***\n", nerrors);
        else{
            printf("===================================\n");
            printf("PHDF5 tests finished with no errors\n");
            printf("===================================\n");
        }
    }
    if (docleanup)
        cleanup();
    MPI_Finalize();

    return(nerrors);
}

#else /* H5_HAVE_PARALLEL */
/* dummy program since H5_HAVE_PARALLE is not configured in */
    int
main(void)
{
    printf("No PHDF5 example because parallel is not configured in\n");
    return(0);
}
#endif /* H5_HAVE_PARALLEL */
