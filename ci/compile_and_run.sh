#!/bin/bash
set -e
cd $1
rm -rf build
mkdir build && cd build
ARGS=(
  -D SCALAPACK_LIBRARY=/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so.2.0
  )
cmake "${ARGS[@]}" ../
make -j8
ctest --no-compress-output -T Test

exit 0
