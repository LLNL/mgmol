#/bin/bash
## An example script to build on ONRL condo systems (CADES).

#Before compiling, load the following modules:
source scripts/modules.condo

MGMOL_ROOT=`pwd`

INSTALL_DIR=${MGMOL_ROOT}/install
mkdir -p ${INSTALL_DIR}

BUILD_DIR=${MGMOL_ROOT}/build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

SCALAPACK_DIR=/home/q8j/Software/ScaLapack/scalapack-2.2.2

# call cmake 
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=mpiCC \
      -DCMAKE_Fortran_COMPILER=mpif77 \
      -DMGMOL_USE_HDF5P=OFF \
      -DMGMOL_WITH_CLANG_FORMAT=ON \
      -DCMAKE_PREFIX_PATH=${HOME}/bin \
      -DSCALAPACK_LIBRARY="${SCALAPACK_DIR}/lib/libscalapack.a;/lib64/libgfortran.so.5" \
      -DMPIEXEC_EXECUTABLE=${OPENMPI_DIR}/bin/mpiexec \
      ..

# call make install
make
make install
