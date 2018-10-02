#Assumes the following packages have been loaded:
#module swap PrgEnv-pgi PrgEnv-gnu
#module load cray-hdf5-parallel
#module load cray-libsci
#module load boost/1.62.0

#platform
PLAT = titan

MPCC    = CC
CXX     = CC

OPT_FLAGS = -O2 -DNDEBUG -g
BLASDIR = $(CRAY_LIBSCI_BASE_DIR)/intel/15.0/x86_64
HDF_DIR = $(HDF5_DIR)
LAPACKLIB  = -lsci_gnu_mp
COMMONLIBS = -lboost_program_options -openmp
CNAME = icpc

NAT_COMPLIB = 

DEBUG_FLAGS = -O0 -g
EXTRACXXFLAGS =

FFLAGS= -O -DNDEBUG

MACHINE_FLAGS = -DADD_ -DHAVE_BOOST

INCLUDE = -I$(HDF_DIR)/include \
	-I$(OLCF_BOOST_ROOT)/include
          
MPI_LIB =

HDF_LIB = $(HDF_DIR)/lib/libhdf5_hl.a $(HDF_DIR)/lib/libhdf5.a -lz

LD_LIBRARY_PATH = -L$(OLCF_BOOST_ROOT)/lib -L$(CRAY_LIBSCI_PREFIX_DIR)/lib

#mixed precision lib
MPBLASLIB =

LIBS= $(LAPACKLIB)

PLIBS = $(HDF_LIB)  

DEFS =
