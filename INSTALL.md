Building MGmol using CMake:
==========================

mkdir build

cd build

cmake ..

make

Third-party libraries will be detected by CMake based on what is in your
environment variables. You can also specify various options to the cmake
command. Examples are provided in the scripts directory.

Third part libraries/packages required:
--------------------------------------
Scalapack

HDF5

Blas

Lapack

Boost

MPI

Code formatting using clang:
---------------------------

To enable code formatting using cmake, use the cmake option
-DMGMOL_WITH_CLANG_FORMAT=ON and specify a path for clang-format 6.0
with -DCMAKE_PREFIX_PATH=/path/to/clang-format.
If you do not have clang-format6.0 already installed, you can download a
precompiled version for Linux at:
https://github.com/dealii/dealii/releases/download/v9.0.0/clang-format-6-linux.tar.gz

Then you can format the code using `make format`
