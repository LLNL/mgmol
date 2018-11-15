#Load following modules:
#module load PE-gnu/3.0
#module load boost
#module load openBLAS
#module load scalapack
#module load hdf5-parallel

#platform
PLAT = cades

##############################################

MPCC    = mpiCC
CXX     = mpiCC

OPT_FLAGS = -O2 -DNDEBUG -g

COMMONLIBS = -lboost_program_options 
CNAME = gnu

NAT_COMPLIB = 

DEBUG_FLAGS = -O0 -g
EXTRACXXFLAGS =

FFLAGS= -O -DNDEBUG

MACHINE_FLAGS = -DADD_ -DHAVE_BOOST

INCLUDE = $(HDF5_INC) \
	  $(LAPACK_INC) $(OPENBLAS_INC) \
	-I$(BOOST_DIR)/include
          
MPI_LIB =

HDF_LIB = $(HDF5_LIB) -lhdf5 -lhdf5_hl
SCALAPACKLIB = $(SCALAPACK_LIB) -lscalapack -lgfortran

LD_LIBRARY_PATH = $(BOOST_LIB)

#mixed precision lib
MPBLASLIB =

LIBS= $(LAPACK_LIB) $(OPENBLAS_LIB) -lopenblas

PLIBS = $(SCALAPACKLIB) $(HDF_LIB)

DEFS =
