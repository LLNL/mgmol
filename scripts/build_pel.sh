##! /bin/csh -f
## An example script to build on LLNL Peloton systems.
## For now, this script assumes intel/ mkl libraries are being used.

#Before compiling, load the following modules:
#module load intel
#module load hdf5-parallel
#module load mkl

# load some modules
module load cmake/3.8.2
module load intel
module load mkl
module load hdf5-parallel

# set some environment variables. Set them explicitly or use loaded module path (preferred)
# Here we use an explicit path for scalapack to be consistent with the path for the blas libraries and avoid 
# benign cmake warnings 
#setenv HDF5_ROOT /usr/tce/packages/hdf5/hdf5-parallel-1.10.2-intel-18.0.1-mvapich2-2.2/
setenv SCALAPACK_ROOT /usr/tce/packages/mkl/mkl-2019.0/
setenv HDF5_ROOT ${HDF5}

# We need to define the cmake blas vendor option here to find the right one.
set BLAS_VENDOR = Intel10_64lp

# manually set the location of BLACS libraries for scalapack
set BLACS_LIB = ${SCALAPACK_ROOT}/lib
#set BLACS_LIB = ${MKLROOT}/lib/intel64/

set MGMOL_ROOT = `pwd`

set INSTALL_DIR = ${MGMOL_ROOT}/mgmol_install
mkdir -p ${INSTALL_DIR}

set BUILD_DIR = ${MGMOL_ROOT}/mgmol_build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake 
#cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_CXX_COMPILER=mpic++ -DBLA_VENDOR=${BLAS_VENDOR} -DSCALAPACK_LIBRARY=${MKLROOT}/lib/intel64/#libmkl_scalapack_lp64.so -DSCALAPACK_INCLUDE_DIR=${MKLROOT}/include/ -DSCALAPACK_BLACS_LIBRARY=${BLACS_LIB}/libmkl_blacs_intelmpi_lp64.so ..
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_CXX_COMPILER=mpic++ -DBLA_VENDOR=${BLAS_VENDOR} -DSCALAPACK_BLACS_LIBRARY=${BLACS_LIB}/libmkl_blacs_intelmpi_lp64.so ..

# call make install
make -j 
make install
