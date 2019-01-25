#platform
PLAT = summit

##############################################
#GNU compiler
##############################################

MPCC    = mpicxx
CXX     = mpicxx
OPT_FLAGS = -O2 -DNDEBUG -g

BLASDIR = $(OLCF_NETLIB_LAPACK_ROOT)
BLASLIBDIR = $(BLASDIR)/lib
HDF_DIR = $(OLCF_HDF5_ROOT)
BOOST_DIR = $(OLCF_BOOST_ROOT)
SCALAPACK_DIR = $(OLCF_NETLIB_SCALAPACK_ROOT)

LAPACKLIB  = -L$(OLCF_NETLIB_LAPACK_ROOT)/lib64 -llapack -lblas
SCALAPACKLIB = -L$(SCALAPACK_DIR)/lib -lscalapack

COMMONLIBS = -lboost_program_options

CNAME = icpc

NAT_COMPLIB = 

DEBUG_FLAGS = -O0 -g -Wuninitialized -ftrapuv
EXTRACXXFLAGS =

FFLAGS= -O -DNDEBUG

MACHINE_FLAGS = -DADD_

INCLUDE = -I$(HDF_DIR)/include \
	-I$(BLASDIR)/include \
	-I$(BOOST_DIR)/include
          
MPI_LIB =

HDF_LIB = -L$(HDF_DIR)/lib -lhdf5 -lhdf5_cpp -lhdf5_hl

MGMOL_LD_LIBRARY_PATH = -L$(BOOST_DIR)/lib -L$(BLASLIBDIR) -Wl,-rpath,$(BLASLIBDIR) \
                  -L$(OLCF_NETLIB_LAPACK_ROOT)/lib64 -Wl,-rpath,$(OLCF_NETLIB_LAPACK_ROOT)/lib64 \
                  -L$(SCALAPACK_DIR)/lib -Wl,-rpath,$(SCALAPACK_DIR)/lib

#mixed precision lib
MPBLASLIB =

LIBS= $(LAPACKLIB) 

PLIBS = $(SCALAPACKLIB) $(HDF_LIB)  

DEFS =
