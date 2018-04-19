#platform
PLAT = bgq

SEQ_ROOT = /bgsys/drivers/ppcfloor/

BLASDIR       = /usr/local/tools/essl/5.1
HDF_DIR       = /usr/local/tools/hdf5/hdf5-1.8.16/parallel_xl
SCALAPACK_DIR = /usr/local/tools/scalapack/lib
BOOST_DIR     = /usr/global/tools/adept/boost-1.47.0/bgqos_0


#compilers
MPCC  = $(SEQ_ROOT)/comm/xl/bin/mpixlcxx_r
MPGCC = /usr/local/bin/mpig++-latestmpi

MACHINE_FLAGS = -DBGQ -DHAVE_BOOST -DUSE_BARRIERS

INCLUDE = -I../include -I$(HDF_DIR)/include \
          -I$(BLASDIR)/include -I$(BOOST_DIR)/include

CFLAGS_BASE =   

DEBUG_FLAGS   = $(CFLAGS_BASE) -g -O0 -qsmp=omp:noopt -qcheck=all -qminimaltoc
OPT_FLAGS     = $(CFLAGS_BASE) -g -DNDEBUG -qstrict -qsmp=omp:auto -O3 -qinfo=cns:uni:ini:pro:ord -qminimaltoc
GCC_OPT_FLAGS = $(CFLAGS_BASE) -g -DNDEBUG -O3
EXTRACXXFLAGS = 

COMMONLIBS = -L/opt/ibmcmp/xlsmp/bg/3.1/bglib64 \
   -L/opt/ibmcmp/xlmass/bg/7.3/bglib64 \
   -L$(BOOST_DIR)/lib/boost-1.47.0 -L/opt/ibmcmp/xlf/bg/14.1/bglib64 \
   -lxlsmp -lxlf90_r -lxlopt -lxl \
   -lboost_program_options-mt -lxlfmath -lm

#compiler specific lib
NAT_COMPLIB     = -lxlomp_ser -qsmp=omp
GCC_COMPLIB = -ldl

MPI_LIB =

HDF_LIB = $(HDF_DIR)/lib/libhdf5_hl.a $(HDF_DIR)/lib/libhdf5.a

ARPACK_LIB=

# Blacs libraries
BLACSdir      = /usr/local/tools/blacs/lib
BLACSFINIT    = $(BLACSdir)/libblacsF77init.a
BLACSCINIT    = $(BLACSdir)/libblacsCinit.a
BLACSLIB      = $(BLACSdir)/libblacs.a

CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)

# Scalapack libraries
SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

LD_LIBRARY_PATH = -L$(BLASDIR)/lib

LAPACKLIB   = /usr/local/tools/lapack/lib/liblapack.a -lesslsmpbg
#LAPACKLIB   = /usr/local/tools/lapack/lib/liblapack.a -lesslbg

#mixed precision lib
MPBLASLIB =
MEMPLIB = -L/usr/local/tools/memP/lib -lmemP -ldl -Wl,-zmuldefs

LIBS= $(LIBDISTMATRIX) $(LAPACKLIB) 

PLIBS = $(SCALAPACKLIB) $(ARPACK_LIB) \
        $(CBLACSLIB) $(LAPACKLIB) $(HDF_LIB)  

DEFS = 
