#!/bin/bash
set -e
cd $1
rm -rf build
mkdir build && cd build
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
