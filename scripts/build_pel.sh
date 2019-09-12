##! /bin/csh -f
## An example script to build on LLNL Peloton systems.
## For now, this script assumes intel/ mkl libraries are being used.

# load some modules
module load cmake
module load intel/19.0.4 
module load mkl
module load hdf5-parallel
module load boost

# set some environment variables. Set them explicitly or use loaded module path (preferred)
# Here we use an explicit path for scalapack to be consistent with the path for the blas libraries and avoid 
# benign cmake warnings 
setenv SCALAPACK_ROOT /usr/tce/packages/mkl/mkl-2019.0/
setenv HDF5_ROOT ${HDF5}

# We need to define the cmake blas vendor option here to find the right one.
set BLAS_VENDOR = Intel10_64lp

# manually set the location of BLACS libraries for scalapack
set BLACS_LIB = ${SCALAPACK_ROOT}/lib

set MGMOL_ROOT = `pwd`

set INSTALL_DIR = ${MGMOL_ROOT}/mgmol_install
mkdir -p ${INSTALL_DIR}

set BUILD_DIR = ${MGMOL_ROOT}/mgmol_build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake 
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_CXX_COMPILER=mpic++ \
      -DMPIEXEC_NUMPROC_FLAG="-n" \
      -DBLA_VENDOR=${BLAS_VENDOR} \
      -DSCALAPACK_BLACS_LIBRARY=${BLACS_LIB}/libmkl_blacs_intelmpi_lp64.so \
      ..

# call make install
make -j 
make install
