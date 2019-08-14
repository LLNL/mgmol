#/bin/bash

export BLAS_VENDOR=OpenBLAS

MGMOL_ROOT=`pwd`

INSTALL_DIR=${MGMOL_ROOT}/mgmol_install
mkdir -p ${INSTALL_DIR}

BUILD_DIR=${MGMOL_ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake 
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_CXX_COMPILER=mpiCC \
      -DBLAS_LIBRARIES=/usr/lib/x86_64-linux-gnu/openblas/libblas.so.3 \
      -DLAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/openblas/liblapack.so.3 \
      -DBLA_VENDOR=OpenBLAS \
      -DCMAKE_Fortran_COMPILER=mpif77 \
      -DSCALAPACK_ROOT=/usr/lib/x86_64-linux-gnu \
      -DMPIEXEC_PREFLAGS="-report-bindings;--map-by;core;-bind-to;core" \
      ..

# call make install
make -j 4
make install
