# $Id:$

#platform
PLAT = pel

##############################################
#Intel compiler
##############################################

MPCC    = /usr/local/bin/mpiicpc-14.0.174
CXX     = /usr/local/bin/icpc-14.0.174
OPT_FLAGS = -O2 -DNDEBUG -Zp16 -g -Wunused-function -Wunused-variable -fp-model precise -fp-model source
BLASDIR = /opt/intel-14.0/mkl
BLASLIBDIR = $(BLASDIR)/lib/intel64
HDF_DIR = /usr/local/tools/hdf5-intel-parallel-1.8.8
HDF_DIR_SERIAL = /usr/local/tools/hdf5-intel-serial-1.8.8
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

NAT_COMPLIB = -lifcore -openmp
GCC_COMPLIB = -Wl,-rpath,/usr/local/tools/ifort-14.0.174/lib -fopenmp

DEBUG_FLAGS = -O0 -g -Wuninitialized -ftrapuv
EXTRACXXFLAGS =

FFLAGS= -O -DNDEBUG

BOOST_DIR      = /usr/local/tools/boost-nompi

LIBXC_DIR      = /usr/gapps/mgmol/Libxc/1.2.0/build

ARPACK_DIR     = $(HOME)/ARPACK

MACHINE_FLAGS = -DADD_ -DHAVE_BOOST -DUSE_CNR

INCLUDE = -I/usr/local/include -I$(HDF_DIR)/include \
	-I$(BLASDIR)/include \
	-I$(BOOST_DIR)/include \
	-I$(LIBXC_DIR)/include
          
MPI_LIB =

HDF_LIB = $(HDF_DIR)/lib/libhdf5_hl.a $(HDF_DIR)/lib/libhdf5.a -lz
HDF_LIB_SERIAL = $(HDF_DIR_SERIAL)/lib/libhdf5.a -lz

ARPACK_LIB=$(ARPACK_DIR)/lib/parpack_$(PLAT).a $(ARPACK_DIR)/lib/libarpack_$(PLAT).a
ARPACK_LIB=

#
#  BLACS setup.  All version need the debug level (0 or 1),
#  and the directory where the BLACS libraries are
#
BLACSDBGLVL   = 0
BLACSdir      = /usr/gapps/mgmol/lib

#
#  MPI setup; uncomment and tailor to your system if using MPIBLACS
#  Will need to comment out the default native BLACS setup below below
#
#USEMPI        = -DUsingMpiBlacs
BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a

#
#  system primitive BLACS setup, comment out if using MPI
#
CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)

#
#  The name of the libraries to be linked to
#
#SCALAPACKLIB  = /usr/gapps/mgmol/lib/scalapack1.8_$(PLAT).a
SCALAPACKLIB  = /opt/intel-14.0/mkl/lib/intel64/libmkl_scalapack_lp64.a \
                /opt/intel-14.0/mkl/lib/intel64/libmkl_blacs_lp64.a 

LD_LIBRARY_PATH = -L$(BOOST_DIR)/lib -L$(BLASLIBDIR) -Wl,-rpath,$(BLASLIBDIR)

#mixed precision lib
MPBLASLIB =

LIBS= $(LIBDISTMATRIX) $(LAPACKLIB) 

PLIBS = $(SCALAPACKLIB) $(ARPACK_LIB) $(CBLACSLIB) $(LAPACKLIB) $(HDF_LIB)  

DEFS =
