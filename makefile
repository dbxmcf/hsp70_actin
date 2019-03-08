PGI_HDF5_BIN=/project/fchen14/packages/hdf5-1.10.5/install/bin

PGI_H5CC=$(PGI_HDF5_BIN)/h5pcc
CFLAGS_PGI_OMP=-mp -fast -Mipa=fast,inline #-Minfo=all
BIN=pomp 

#PCB_MPI_SRC=src/c/ph5ex.c
PCB_MPI_SRC=src/c/ph5_cmb_blks.c

pomp: $(PCB_MPI_SRC)
	$(PGI_H5CC) $(CFLAGS_PGI_OMP) $(PCB_MPI_SRC) -o pgi.mpi.$@.out

all: $(BIN)

clean:
	rm -rf pgi.mpi.omp.out 
