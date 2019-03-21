#/bin/bash
## An example script to build on ONRL condo systems (CADES).
## For now, this script assumes intel/ mkl libraries are being used.

#Before compiling, load the following modules:
source scripts/modules.condo

# set some environment variables. Set them explicitly or use loaded module path (preferred)
# Here we use an explicit path for scalapack to be consistent with the path for the blas libraries and avoid 
export SCALAPACK_ROOT=${MKLROOT}
export HDF5_ROOT=${HDF5_PARALLEL_DIR}

# We need to define the cmake blas vendor option here to find the right one.
BLAS_VENDOR=Intel10_64lp

# manually set the location of BLACS libraries for scalapack
#set BLACS_LIB = ${SCALAPACK_ROOT}/lib
BLACS_LIB=${MKLROOT}/lib/intel64

MGMOL_ROOT=`pwd`

INSTALL_DIR=${MGMOL_ROOT}/mgmol_install
mkdir -p ${INSTALL_DIR}

BUILD_DIR=${MGMOL_ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake 
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_CXX_COMPILER=mpiCC -DBLA_VENDOR=${BLAS_VENDOR} -DSCALAPACK_BLACS_LIBRARY=${BLACS_LIB}/libmkl_blacs_openmpi_lp64.so ..

# call make install
make -j 
make install
