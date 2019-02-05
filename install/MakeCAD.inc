#Load following modules:
#module load PE-intel/3.0
#module load hdf5_parallel/1.8.21
#module load boost/1.67.0-pe3
#module load mkl

#platform
PLAT = cades

##############################################

MPCC    = mpiCC
CXX     = mpiCC

OPT_FLAGS = -O2 -DNDEBUG -g 

COMMONLIBS = -lboost_program_options 
CNAME = 

NAT_COMPLIB = 

DEBUG_FLAGS = -O0 -g
EXTRACXXFLAGS =

FFLAGS= -O -DNDEBUG

MACHINE_FLAGS = -DADD_

INCLUDE = -I $(HDF5_PARALLEL_DIR)/include \
	-I$(BOOST_DIR)/include
          
MPI_LIB =

HDF_LIB = -L $(HDF5_PARALLEL_DIR)/lib -lhdf5 -lhdf5_hl
SCALAPACKLIB = -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64

MGMOL_LD_LIBRARY_PATH = 

#mixed precision lib
MPBLASLIB =

LIBS= $(MKL_LIB) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
PLIBS = $(SCALAPACKLIB) $(HDF_LIB)

DEFS =
