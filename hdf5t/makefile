INTEL_HDF5_BIN=/usr/local/packages/hdf5/1.8.12/INTEL-140-MVAPICH2-2.0/bin
PGI_HDF5_BIN=/project/fchen14/packages/hdf5-1.10.4/install/bin

INTEL_H5CC=$(INTEL_HDF5_BIN)/h5pcc
PGI_H5CC=$(PGI_HDF5_BIN)/h5pcc
#CFLAGS=-acc -Minfo=accel -ta=nvidia,time $(INC)
CFLAGS_INTEL_OMP=-openmp
CFLAGS_PGI_OMP=-mp -fast -Minfo=all
CFLAGS_INTEL_OMP=-mp -fast -Minfo=all
CFLAGS_ACC=-acc -fast -Minfo=all -ta=nvidia,time

PCB_MPI_OMP=ph5_cmb_blks.c

pgi.mpi.omp: $(PCB_MPI_OMP)
	$(PGI_H5CC) $(CFLAGS_PGI_OMP) $(PCB_MPI_OMP) -o $@.out

intel.mpi.omp: $(PCB_MPI_OMP)
	$(INTEL_H5CC) $(CFLAGS_INTEL_OMP) $(PCB_MPI_OMP) -o $@.out

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
	rm -rf pgi.omp.out intel.omp.out
