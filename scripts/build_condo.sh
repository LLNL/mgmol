#/bin/bash
## An example script to build on ONRL condo systems (CADES).
## This script assumes intel/ mkl libraries are being used.

#Before compiling, load the following modules:
source scripts/modules.condo

# set some environment variables using loaded module path
export SCALAPACK_ROOT=${MKLROOT}

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
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=mpiCC \
      -DCMAKE_Fortran_COMPILER=mpif77 \
      -DHDF5_LIBRARIES=${HDF5_DIR}/lib/libhdf5.so \
      -DHDF5_HL_LIBRARIES=${HDF5_DIR}/lib/libhdf5_hl.so \
      -DHDF5_INCLUDE_DIRS=${HDF5_DIR}/include \
      -DBLA_VENDOR=${BLAS_VENDOR} \
      -DMGMOL_WITH_CLANG_FORMAT=ON \
      -DCMAKE_PREFIX_PATH=${HOME}/bin \
      -DMPIEXEC_EXECUTABLE=${OPENMPI_DIR}/bin/mpiexec \
      -DSCALAPACK_BLACS_LIBRARY=${BLACS_LIB}/libmkl_blacs_openmpi_lp64.so \
      ..

# call make install
make
make install
