#/bin/bash
## An example script to build on ONRL condo systems (CADES).
## This script assumes intel/ mkl libraries are being used.

#Before compiling, load the following modules:
source scripts/modules.condo

# set some environment variables using loaded module path
export SCALAPACK_ROOT=${MKLROOT}
export HDF5_ROOT=${HDF5_PARALLEL_DIR}

# We need to define the cmake blas vendor option here to find the right one.
BLAS_VENDOR=Intel10_64lp

# manually set the location of BLACS libraries for scalapack
BLACS_LIB=${MKLROOT}/lib/intel64

MGMOL_ROOT=`pwd`

INSTALL_DIR=${MGMOL_ROOT}/mgmol_install
mkdir -p ${INSTALL_DIR}

BUILD_DIR=${MGMOL_ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake 
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_CXX_COMPILER=mpiCC \
      -DCMAKE_Fortran_COMPILER=mpif77 \
      -DBLA_VENDOR=${BLAS_VENDOR} \
      -DMPIEXEC_EXECUTABLE=${OPENMPI_DIR}/bin/mpiexec \
      -DSCALAPACK_BLACS_LIBRARY=${BLACS_LIB}/libmkl_blacs_openmpi_lp64.so \
      ..

# call make install
make -j 
make install
