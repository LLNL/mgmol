#/bin/bash

#Before compiling, load the following modules:
source scripts/modules.condo-mod

# We need to define the cmake blas vendor option here to find the right one.
BLAS_VENDOR=OpenBLAS

MGMOL_ROOT=`pwd`

INSTALL_DIR=${MGMOL_ROOT}/install
mkdir -p ${INSTALL_DIR}

BUILD_DIR=${MGMOL_ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# call cmake 
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=mpiCC \
      -DCMAKE_Fortran_COMPILER=mpif77 \
      -DBLA_VENDOR=${BLAS_VENDOR} \
      -DMGMOL_USE_HDF5P=OFF \
      -DMGMOL_WITH_CLANG_FORMAT=ON \
      -DCMAKE_PREFIX_PATH=${HOME}/bin \
      -DSCALAPACK_LIBRARY="${SCALAPACK_DIR}/lib/libscalapack.a;/lib64/libgfortran.so.3" \
      -DMPIEXEC_EXECUTABLE=${OPENMPI_DIR}/bin/mpiexec \
      ..

# call make install
make
make install
