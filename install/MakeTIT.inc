#Assumes the following packages have been loaded:
#module swap PrgEnv-pgi PrgEnv-gnu
#module load cray-hdf5-parallel
#module load cray-libsci
#module load boost/1.62.0
#module load scorep

#platform
PLAT = titan

#PREF = scorep

MPCC    = $(PREF) CC
CXX     = $(PREF) CC

OPT_FLAGS = -O2 -DNDEBUG -g -fopenmp -ftree-vectorize

BLASDIR = $(CRAY_LIBSCI_BASE_DIR)/intel/15.0/x86_64
HDF_DIR = $(HDF5_DIR)
LAPACKLIB  = -lsci_gnu_mp

COMMONLIBS = -lboost_program_options -fopenmp

CNAME = icpc

NAT_COMPLIB = 

DEBUG_FLAGS = -O0 -g
EXTRACXXFLAGS =

FFLAGS= -O -DNDEBUG

MACHINE_FLAGS = -DADD_

INCLUDE = -I$(HDF_DIR)/include \
	-I$(OLCF_BOOST_ROOT)/include
          
MPI_LIB =

HDF_LIB = $(HDF_DIR)/lib/libhdf5_hl.a $(HDF_DIR)/lib/libhdf5.a -lz

MGMOL_LD_LIBRARY_PATH = -L$(OLCF_BOOST_ROOT)/lib -L$(CRAY_LIBSCI_PREFIX_DIR)/lib

#mixed precision lib
MPBLASLIB =

LIBS= $(LAPACKLIB)

PLIBS = $(HDF_LIB)  

DEFS =
