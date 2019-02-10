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
#ifndef PH5_CMB_BLKS
#define PH5_CMB_BLKS
#include <assert.h>
#include "hdf5.h"
#include <string.h>
#include <stdlib.h>

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

/* global variables */
int nerrors = 0;				/* errors count */
#ifndef PATH_MAX
#define PATH_MAX    512
#endif  /* !PATH_MAX */
char    testfiles[2][PATH_MAX];


int mpi_size, mpi_rank;				/* mpi variables */

/* option flags */
int verbose = 0;			/* verbose, default as no. */
int doread=1;				/* read test */
int dowrite=1;				/* write test */
int docleanup=0;			/* cleanup */

/* Prototypes */
void slab_set(hsize_t start[], hsize_t count[], hsize_t stride[], int mode);
void dataset_fill(hsize_t start[], hsize_t count[], hsize_t stride[], DATATYPE * dataset);
void dataset_print(hsize_t start[], hsize_t count[], hsize_t stride[], DATATYPE * dataset);
int dataset_vrfy(hsize_t start[], hsize_t count[], hsize_t stride[], DATATYPE *dataset, DATATYPE *original);
void phdf5writeInd(char *filename);
void phdf5readInd(char *filename);
void phdf5writeAll(char *filename);
void phdf5readAll(char *filename);
void test_split_comm_access(char filenames[][PATH_MAX]);
int parse_options(int argc, char **argv);
void usage(void);
int mkfilenames(char *prefix);
void cleanup(void);
#endif