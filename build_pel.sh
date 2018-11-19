##! /bin/csh -f
## An example script to build on LLNL Peloton systems.
## For now, this script assumes intel/ mkl libraries are being used.

# load some modules
module load cmake/3.8.2
module load intel
module load mkl
module load hdf5-parallel

# set some environment variables
setenv HDF5_ROOT /usr/tce/packages/hdf5/hdf5-parallel-1.10.2-intel-18.0.1-mvapich2-2.2/
setenv SCALAPACK_ROOT /usr/tce/packages/mkl/mkl-2019.0/

# We need to define the cmake blas vendor option here, to find the right one.
set BLAS_VENDOR = Intel10_64lp

# manually set the location of BLACS libraries for scalapack
set BLACS_LIB = ${SCALAPACK_ROOT}/lib

set INSTALL_DIR = mgmol_install
mkdir -p ${INSTALL_DIR}

set BUILD_DIR = mgmol_build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake 
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_CXX_COMPILER=mpic++ -DSCALAPACK_BLACS_LIBRARIES=${BLACS_LIB}/libmkl_blacs_intelmpi_lp64.so -DBLA_VENDOR=${BLAS_VENDOR} ..

# call make install
make -j install 
