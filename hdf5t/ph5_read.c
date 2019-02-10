#include "ph5_cmb_blks.h"

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
    if ((SPACE1_DIM1 % mpi_size) || (SPACE1_DIM2 % mpi_size)){
	printf("DIM1(%d) and DIM2(%d) must be multiples of processes (%d)\n",
	    SPACE1_DIM1, SPACE1_DIM2, mpi_size);
	nerrors++;
	goto finish;
    }

    if (parse_options(argc, argv) != 0)
	goto finish;

    /* show test file names */
    if (mpi_rank == 0){
	n = sizeof(testfiles)/sizeof(testfiles[0]);
	printf("Parallel test files are:\n");
	for (i=0; i<n; i++){
	    printf("   %s\n", testfiles[i]);
	}
    }

    if (doread){
	MPI_BANNER("testing PHDF5 dataset independent read...");
	phdf5readInd(testfiles[0]);
	MPI_BANNER("testing PHDF5 dataset collective read...");
	phdf5readAll(testfiles[1]);
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