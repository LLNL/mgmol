#/bin/bash

#Before compiling, load the following modules:
source scripts/modules.condo-llvm

MGMOL_ROOT=`pwd`

INSTALL_DIR=${MGMOL_ROOT}/install
mkdir -p ${INSTALL_DIR}

BUILD_DIR=${MGMOL_ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

SCALAPACK_DIR=$OLCF_NETLIB_SCALAPACK_ROOT

# call cmake 
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=mpiCC \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_Fortran_COMPILER=mpif77 \
      -DMGMOL_USE_HDF5P=OFF \
      -DMGMOL_WITH_CLANG_FORMAT=OFF \
      -DMPIEXEC_EXECUTABLE=${OPENMPI_DIR}/bin/mpiexec \
      ..

# call make install
make -j4
make install
