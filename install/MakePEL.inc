#Before compiling, load the following modules:
#module load intel
#module load hdf5-parallel
#module load mkl

#platform
PLAT = pel

##############################################
#Intel compiler
##############################################

MPCC    = mpic++
CXX     = mpic++
OPT_FLAGS = -O2 -DNDEBUG -Zp16 -g -Wunused-function -Wunused-variable -fp-model precise -fp-model source
BLASDIR = $(MKLROOT)
BLASLIBDIR = $(BLASDIR)/lib/intel64
HDF_DIR = $(HDF5)
LAPACKLIB  = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
COMMONLIBS = -lboost_program_options
CNAME = icpc

##############################################
#gnu compiler
##############################################

MPGCC = /usr/local/bin/mpig++
GCC_OPT_FLAGS   =  -O2 -DNDEBUG -g
#BLASDIR = /opt/intel-14.0/mkl
#BLASLIBDIR = $(BLASDIR)/lib/intel64
##BLASDIR = /usr/lib64/atlas
##BLASLIBDIR = /usr/lib64/atlas/lib
#HDF_DIR = /usr/local/tools/hdf5-gnu-parallel-1.8.8
#HDF_DIR_SERIAL = /usr/local/tools/hdf5-gnu-serial-1.8.8
#LAPACKLIB  = -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core
##LAPACKLIB  = -llapack -L/usr/local/tools/ifort-14.0.174/lib -lifcore -limf -lirc -lsvml
#CNAME = g++

##############################################

NAT_COMPLIB = -lifcore -qopenmp
GCC_COMPLIB = -Wl,-rpath,/usr/local/tools/ifort-14.0.174/lib -fopenmp

DEBUG_FLAGS = -O0 -g -Wuninitialized -ftrapuv
EXTRACXXFLAGS =

FFLAGS= -O -DNDEBUG

#BOOST_DIR      = /usr/local/tools/boost-nompi

LIBXC_DIR      = /usr/gapps/mgmol/Libxc/1.2.0/build

ARPACK_DIR     = 

MACHINE_FLAGS = -DADD_ -DUSE_CNR

INCLUDE = -I$(HDF_DIR)/include \
	-I$(BLASDIR)/include
          
MPI_LIB =

HDF_LIB = $(HDF_DIR)/lib/libhdf5_hl.a $(HDF_DIR)/lib/libhdf5.a -lz
HDF_LIB_SERIAL = $(HDF_DIR_SERIAL)/lib/libhdf5.a -lz

ARPACK_LIB=

#
#  The name of the libraries to be linked to
#
SCALAPACKLIB  = $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a \
                $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a 

LD_LIBRARY_PATH = -L$(BLASLIBDIR) -Wl,-rpath,$(BLASLIBDIR)

#mixed precision lib
MPBLASLIB =

LIBS= $(LAPACKLIB) 

PLIBS = $(SCALAPACKLIB) $(HDF_LIB)

DEFS =
