##! /bin/csh -f
## An example script to build on LLNL Peloton systems.
## For now, this script assumes intel/ mkl libraries are being used.

# load some modules
source scripts/modules.quartz

# set some environment variables. Set them explicitly or use loaded module path (preferred)
# Here we use an explicit path for scalapack to be consistent with the path for the blas libraries and avoid
# benign cmake warnings
#setenv SCALAPACK_ROOT /usr/tce/packages/mkl/mkl-2020.0/lib
setenv SCALAPACK_ROOT ${MKLROOT}
setenv HDF5_ROOT /usr/tce/packages/hdf5/hdf5-1.14.0-mvapich2-2.3.6-intel-2022.1.0

# We need to define the cmake blas vendor option here to find the right one.
set BLAS_VENDOR = Intel10_64lp

# manually set the location of BLACS libraries for scalapack
set BLACS_LIB = ${SCALAPACK_ROOT}/lib/intel64

set MGMOL_ROOT = `pwd`

set INSTALL_DIR = ${MGMOL_ROOT}/install_quartz
mkdir -p ${INSTALL_DIR}

set BUILD_DIR = ${MGMOL_ROOT}/build_quartz
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# clone the libROM GitHub repo in BUILD_DIR
set USE_LIBROM="On"
set LIBROM_PATH = ${BUILD_DIR}/libROM
git clone https://github.com/LLNL/libROM

# call cmake
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_CXX_COMPILER=mpic++ \
      -DCMAKE_Fortran_COMPILER=mpif77 \
      -DMPIEXEC_NUMPROC_FLAG="-n" \
      -DBLA_VENDOR=${BLAS_VENDOR} \
      -DSCALAPACK_BLACS_LIBRARY=${BLACS_LIB}/libmkl_blacs_intelmpi_lp64.so \
      -DUSE_LIBROM=${USE_LIBROM} \
      -DLIBROM_PATH=${LIBROM_PATH} \
      -DCMAKE_BUILD_TYPE=DEBUG \
      ..

# call make install
make -j 
make install

#      -DBLAS_LIBRARIES=/usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/lib \
#      -DLAPACK_LIBRARIES=/usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/lib \  