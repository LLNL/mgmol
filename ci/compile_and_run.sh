#!/bin/bash
set -e
cd $1
rm -rf build
mkdir build && cd build
ARGS=(
  -D SCALAPACK_ROOT=/usr/lib/x86_64-linux-gnu
  -D MGMOL_WITH_CLANG_FORMAT=ON
  )
cmake "${ARGS[@]}" ../
make -j8
ctest --no-compress-output -T Test
make format && git diff --exit-code

exit 0
