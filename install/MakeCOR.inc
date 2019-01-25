#platform
PLAT = cori

##############################################
#module swicth PrgEnv-intel/6.0.4 PrgEnv-gnu
#module load gcc/7.1.0
#module load llvm
#module load cray-hdf5-parallel
#module load boost/1.67.0
#module load cray-libsci

##############################################

MPCC    = CC
CXX     = CC
OPT_FLAGS = -O2 -DNDEBUG -g
BLASDIR = $(CRAY_LIBSCI_DIR)/gnu/7.1/haswell
BLASLIBDIR = $BLASDIR/lib
HDF_DIR = $(HDF5_ROOT)
LAPACKLIB  = -lsci_gnu_mp
COMMONLIBS = -lboost_program_options -openmp
CNAME = icpc

DEBUG_FLAGS = -O0 -g -Wuninitialized -ftrapuv
EXTRACXXFLAGS =

FFLAGS= -O -DNDEBUG

MACHINE_FLAGS = -DADD_

INCLUDE = -I$(HDF_DIR)/include \
	-I$(BLASDIR)/include \
	-I$(BOOST_DIR)/include
          
MPI_LIB =

HDF_LIB = $(HDF_DIR)/lib/libhdf5_hl.a $(HDF_DIR)/lib/libhdf5.a -lz

MGMOL_LD_LIBRARY_PATH = -L$(BOOST_DIR)/lib -L$(BLASLIBDIR) -Wl,-rpath,$(BLASLIBDIR)

#mixed precision lib
MPBLASLIB =

LIBS= $(LIBDISTMATRIX) $(LAPACKLIB) 

PLIBS = $(LAPACKLIB) $(HDF_LIB)  

DEFS =
