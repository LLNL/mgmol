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
      -DCMAKE_BUILD_TYPE="Release" \
      -DCMAKE_CXX_COMPILER=mpiCC.openmpi \
      -DBLAS_LIBRARIES=/usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblas.so \
      -DLAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so \
      -DBLA_VENDOR=OpenBLAS \
      -DCMAKE_Fortran_COMPILER=mpif77.openmpi \
      -DSCALAPACK_LIBRARY=/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so.2.1 \
      -DMPIEXEC_EXECUTABLE=/usr/bin/mpirun \
      -DMPIEXEC_NUMPROC_FLAG="-np" \
      -DMPIEXEC_PREFLAGS="--oversubscribe" \
      -DMGMOL_WITH_CLANG_FORMAT=ON \
      -DCMAKE_PREFIX_PATH=${HOME}/bin \
      -D CMAKE_CXX_FLAGS="-Wall -pedantic -Wextra -Wno-cast-function-type" \
      ..

# call make install
make -j 4
make install
