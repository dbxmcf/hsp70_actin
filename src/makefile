INTEL_HDF5_BIN=/usr/local/packages/hdf5/1.8.12/INTEL-140-MVAPICH2-2.0/bin
PGI_HDF5_BIN=/project/fchen14/packages/hdf5-1.10.4/install/bin

INTEL_H5CC=$(INTEL_HDF5_BIN)/h5pcc
PGI_H5CC=$(PGI_HDF5_BIN)/h5pcc
#CFLAGS=-acc -Minfo=accel -ta=nvidia,time $(INC)
CFLAGS_INTEL_OMP=-openmp -O3
CFLAGS_PGI_OMP=-mp -fast -Mipa=fast,inline #-Minfo=all
#CFLAGS_PGI_ACC=-acc -fast -Mipa=fast,inline -Minfo=all -ta=nvidia,time
CFLAGS_PGI_ACC=-acc -fast -Mipa=fast,inline -Mcuda -ta=nvidia:cc3+,cuda9.0,fastmath
#BIN=pgi.mpi.omp intel.mpi.omp
BIN=pomp pacc

PCB_MPI_SRC=ph5_cmb_blks.c

pomp: $(PCB_MPI_SRC)
	$(PGI_H5CC) $(CFLAGS_PGI_OMP) $(PCB_MPI_SRC) -o pgi.mpi.$@.out

pacc: $(PCB_MPI_SRC)
	$(PGI_H5CC) $(CFLAGS_PGI_ACC) $(PCB_MPI_SRC) -o pgi.mpi.$@.out

intel: $(PCB_MPI_SRC)
	$(INTEL_H5CC) $(CFLAGS_INTEL_OMP) $(PCB_MPI_SRC) -o $@.mpi.omp.out

#lap_acc_v0_dr: $(LAP_ACC)
#	$(CC) $(CFLAGS_ACC) $(LAP_ACC_V0_DR) -o $@.out
#
#lap_acc_v1: $(LAP_ACC)
#	$(CC) $(CFLAGS_ACC) $(LAP_ACC_V1) -o $@.out
#
#lap_omp: $(LAP_ACC)
#	$(CC) $(CFLAGS_OMP) $(LAP_ACC_V0) -o $@.out
#
#lap_serial: $(LAP_ACC)
#	$(CC) $(CFLAGS_SER) $(LAP_ACC_V0) -o $@.out

all: $(BIN)

clean:
	rm -rf pgi.mpi.omp.out intel.mpi.omp.out
