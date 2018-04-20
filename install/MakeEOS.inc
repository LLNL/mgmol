# $Id:$

#platform
PLAT = eos

##############################################
#Intel compiler
##############################################

MPCC    = CC
CXX     = CC
OPT_FLAGS = -O2 -DNDEBUG -g
BLASDIR = $(CRAY_LIBSCI_BASE_DIR)/intel/15.0/x86_64
BLASLIBDIR = $(BLASDIR)/lib
HDF_DIR = $(HDF5_DIR)
LAPACKLIB  = -lsci_intel_mp
COMMONLIBS = -lboost_program_options 
CNAME = icpc

NAT_COMPLIB = -lifcore -qopenmp

DEBUG_FLAGS = -O0 -g -Wuninitialized -ftrapuv
EXTRACXXFLAGS =

FFLAGS= -O -DNDEBUG

MACHINE_FLAGS = -DADD_ -DHAVE_BOOST

INCLUDE = -I$(HDF_DIR)/include \
	-I$(BLASDIR)/include \
	-I$(BOOST_DIR)/include
          
MPI_LIB =

HDF_LIB = $(HDF_DIR)/lib/libhdf5_hl.a $(HDF_DIR)/lib/libhdf5.a -lz

LD_LIBRARY_PATH = -L$(BOOST_DIR)/lib -L$(BLASLIBDIR) -Wl,-rpath,$(BLASLIBDIR)

#mixed precision lib
MPBLASLIB =

LIBS= $(LIBDISTMATRIX) $(LAPACKLIB) 

PLIBS = $(LAPACKLIB) $(HDF_LIB)  

DEFS =
