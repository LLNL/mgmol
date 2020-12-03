#! /bin/csh -f

# load some modules
./scripts/modules.summit

# set some environment variables. Set them explicitly or use loaded module path (preferred)

# hdf5
setenv HDF5_ROOT ${OLCF_HDF5_ROOT}

setenv SCALAPACK_ROOT ${OLCF_NETLIB_SCALAPACK_ROOT}

# We need to define the cmake blas vendor option here to find the right one.
set BLA_VENDOR = OpenBLAS

set MGMOL_ROOT = `pwd`

set INSTALL_DIR = ${MGMOL_ROOT}/mgmol_install
mkdir -p ${INSTALL_DIR}

set BUILD_DIR = ${MGMOL_ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

cmake \
      -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_CXX_COMPILER=mpicxx \
      -DMPIEXEC_EXECUTABLE="/sw/summit/xalt/1.2.1/bin/jsrun" \
      -DMPIEXEC_NUMPROCS_FLAG="-n" \
      -DMPIEXEC_PREFLAGS="-a1;-c4;-bpacked:2;-g1" \
      -DBLA_VENDOR=${BLA_VENDOR} \
      -DMGMOL_WITH_COVERAGE=OFF \
      -DMGMOL_WITH_MAGMA=ON \
      -DMGMOL_WITH_OPENMP_OFFLOAD=ON \
      -DCMAKE_CXX_FLAGS="-fopenmp -foffload=nvptx-none -foffload="-O3" -fno-stack-protector" \
      -DCMAKE_PREFIX_PATH="/autofs/nccs-svm1_proj/csc304/magma-2.5.3" \
      -DLAPACK_LIBRARIES=${OLCF_NETLIB_LAPACK_ROOT}/lib64/liblapack.so \
      ..

# call make install
make -j 16
make install
