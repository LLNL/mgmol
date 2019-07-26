##! /bin/csh -f

# load some modules
source scripts/modules.summit
#module load cmake/3.8.2
#module load xl/2018.11.26
#module load essl
#module load lapack/3.8.0-gcc-4.9.3
#module load hdf5-parallel

# set some environment variables. Set them explicitly or use loaded module path (preferred)
#
# set SCALAPACK_ROOT to automatically attempt to find scalapack libraries and headers. This is 
# useful if the libraries/ headers have well-defined names like scalapack.a(.so) and scalapack.h.
# Otherwise provide paths for SCALAPACK_LIBRARY and DSCALAPACK_INCLUDE_DIR (as below). 
#setenv SCALAPACK_ROOT ${LAPACK_DIR}

# hdf5
setenv HDF5_ROOT ${OLCF_HDF5_ROOT}

# We need to define the cmake blas vendor option here to find the right one.
set BLAS_VENDOR = IBMESSL

# manually set the location of BLACS libraries for scalapack
#set BLACS_LIB = ${SCALAPACK_ROOT}/lib
#set BLACS_LIB = ${MKLROOT}/lib/intel64/

set MGMOL_ROOT = `pwd`

set INSTALL_DIR = ${MGMOL_ROOT}/mgmol_install
mkdir -p ${INSTALL_DIR}

set BUILD_DIR = ${MGMOL_ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake 
#      -DSCALAPACK_LIBRARY=/usr/tcetmp/packages/scalapack/scalapack-2.0.2-gfortran-4.8.5/lib/libscalapack.a \
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_CXX_COMPILER=mpicxx \
      -DMPIEXEC_EXECUTABLE="/sw/summit/xalt/1.1.3/bin/jsrun" \
      -DMPIEXEC_NUMPROCS_FLAG="-n" \
      -DMPIEXEC_PREFLAGS="-a1;-c7;-bpacked:2;-g1" \
      -DBLA_VENDOR=${BLAS_VENDOR} \
      -DLAPACK_LIBRARIES=${OLCF_ESSL_ROOT}/lib64/libesslsmp.so \
      -DSCALAPACK_LIBRARY=${OLCF_NETLIB_SCALAPACK_ROOT}/lib/libscalapack.so \
#      -DSCALAPACK_INCLUDE_DIR=${LAPACK_DIR} \
      ..

# call make install
#make
make -j
make install
