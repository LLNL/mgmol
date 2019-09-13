#!/bin/bash
set -e
cd $1
rm -rf build
mkdir build && cd build
# Get rid of Read -1, expected <someNumber>, errno =1 error
# See https://github.com/open-mpi/ompi/issues/4948
export OMPI_MCA_btl_vader_single_copy_mechanism=none

ARGS=(
  -D SCALAPACK_ROOT=/usr/lib/x86_64-linux-gnu
  -D MGMOL_WITH_CLANG_FORMAT=ON
  -D MGMOL_WITH_COVERAGE=ON
  )
cmake "${ARGS[@]}" ../
make -j8
ctest --no-compress-output -T Test
make format && git diff --exit-code

# Code coverage
make coverage
curl -s https://codecov.io/bash -o codecov_bash_uploader
chmod +x codecov_bash_uploader
./codecov_bash_uploader -Z -X gcov -f lcov.info

exit 0
