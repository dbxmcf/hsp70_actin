PGI_HDF5_BIN=/project/fchen14/packages/hdf5-1.10.4/install/bin

PGI_H5CC=$(PGI_HDF5_BIN)/h5pcc
CFLAGS_PGI_OMP=-mp -fast -Mipa=fast,inline #-Minfo=all
CFLAGS_PGI_ACC=-acc -fast -Mipa=fast,inline -Mcuda -ta=nvidia:cc3+,cuda9.0,fastmath
BIN=pomp pacc

#PCB_MPI_SRC=src/c/ph5ex.c
PCB_MPI_SRC=src/c/ph5_cmb_blks.c

pomp: $(PCB_MPI_SRC)
	$(PGI_H5CC) $(CFLAGS_PGI_OMP) $(PCB_MPI_SRC) -o pgi.mpi.$@.out

pacc: $(PCB_MPI_SRC)
	$(PGI_H5CC) $(CFLAGS_PGI_ACC) $(PCB_MPI_SRC) -o pgi.mpi.$@.out

all: $(BIN)

clean:
	rm -rf pgi.mpi.omp.out pgi.mpi.acc.out
