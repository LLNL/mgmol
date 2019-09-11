##! /bin/csh -f

# load some modules
source scripts/modules.summit

# set some environment variables. Set them explicitly or use loaded module path (preferred)

# hdf5
setenv HDF5_ROOT ${OLCF_HDF5_ROOT}

setenv SCALAPACK_ROOT ${OLCF_NETLIB_SCALAPACK_ROOT}

# We need to define the cmake blas vendor option here to find the right one.
set BLA_VENDOR = IBMESSL

set MGMOL_ROOT = `pwd`

set INSTALL_DIR = ${MGMOL_ROOT}/mgmol_install
mkdir -p ${INSTALL_DIR}

set BUILD_DIR = ${MGMOL_ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_CXX_COMPILER=mpicxx \
      -DMPIEXEC_EXECUTABLE="/sw/summit/xalt/1.1.3/bin/jsrun" \
      -DMPIEXEC_NUMPROCS_FLAG="-n" \
      -DMPIEXEC_PREFLAGS="-a1;-c7;-bpacked:2;-g1" \
      -DBLA_VENDOR=${BLA_VENDOR} \
      -DLAPACK_LIBRARIES=${OLCF_ESSL_ROOT}/lib64/libesslsmp.so \
      ..

# call make install
make -j
make install
