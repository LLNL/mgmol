##! /bin/csh -f
## An example script to build on LLNL Peloton systems.
## For now, this script assumes intel/ mkl libraries are being used.
 
# load some modules
source scripts/modules.pel
 
# set some environment variables. Set them explicitly or use loaded module path (preferred)
# Here we use an explicit path for scalapack to be consistent with the path for the blas libraries and avoid 
# benign cmake warnings 
setenv SCALAPACK_ROOT ${MKLROOT}
setenv HDF5_ROOT ${HDF5}
 
# We need to define the cmake blas vendor option here to find the right one.
set BLAS_VENDOR = Intel10_64lp
 
# manually set the location of BLACS libraries for scalapack
set BLACS_LIB = ${SCALAPACK_ROOT}/lib/intel64
 
set MGMOL_ROOT = `pwd`
 
set INSTALL_DIR = ${MGMOL_ROOT}/mgmol_install
mkdir -p ${INSTALL_DIR}
 
set BUILD_DIR = ${MGMOL_ROOT}/mgmol_build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}
 
# call cmake 
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_CXX_COMPILER=mpic++ \
      -DCMAKE_Fortran_COMPILER=mpif77 \
      -DMPIEXEC_NUMPROC_FLAG="-n" \
      -DBLA_VENDOR=${BLAS_VENDOR} \
      -DSCALAPACK_BLACS_LIBRARY=${BLACS_LIB}/libmkl_blacs_intelmpi_lp64.so \
      -DCMAKE_BUILD_TYPE=Release \
      -DMPIEXEC_EXECUTABLE=/usr/bin/srun \
      ..
 
# call make install
make -j 
make install


