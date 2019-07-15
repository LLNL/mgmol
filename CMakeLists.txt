cmake_minimum_required (VERSION 3.8)

project (MGmol C CXX Fortran)

set (CMAKE_CXX_STANDARD 11)

# Specify the location of additional CMAKE modules
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

# Version number
set (PROJECT_VERSION_MAJOR 1)
set (PROJECT_VERSION_MINOR 0)
set(MGMOL_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")

# Prohibit in-source build
IF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "In-source build prohibited.")
ENDIF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

# Use RPATH in install, many mpicc scripts use RPATH so default
# behavior of CMAKE to remove RPATH from installed executables is not
# so good.
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Cmake build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "No user specified build type, using default...")
endif()
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)

# set some default flag options
# GNU
set(GNU_CXX_FLAGS_DEBUG -O0 -g -Wuninitialized)
set(GNU_CXX_FLAGS_RELEASE -g -DNDEBUG -O3)
# INTEL
set(Intel_CXX_FLAGS_DEBUG -O0 -g)
set(Intel_CXX_FLAGS_RELEASE -O2 -DNDEBUG -Zp16 -g -Wunused-function -Wunused-variable -fp-model precise -fp-model source)
#set(Intel_CXX_FLAGS_RELEASE -g -DNDEBUG -O3)
# IBM
set(XL_CXX_FLAGS_DEBUG -g -O0 -qsmp=omp:noopt -qcheck=all -qminimaltoc)
set(XL_CXX_FLAGS_RELEASE -g -DNDEBUG -qstrict -qsmp=omp:auto -O3 -qinfo=cns:uni:ini:pro:ord -qminimaltoc)
# Clang
set(Clang_CXX_FLAGS_DEBUG -g -O0 -Wno-deprecated-register)
set(Clang_CXX_FLAGS_RELEASE -g -DNDEBUG -O3 -qstrict -Wno-deprecated-register)
# default
set(DEFAULT_CXX_FLAGS_DEBUG -g -O0)
set(DEFAULT_CXX_FLAGS_RELEASE -g -DNDEBUG -O3)

if(CMAKE_BUILD_TYPE_UPPER)
  if(CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG" 
    OR CMAKE_BUILD_TYPE_UPPER STREQUAL "RELEASE")
    if(CMAKE_CXX_FLAGS STREQUAL "")
      if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
          OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel"
          OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang"
          OR CMAKE_CXX_COMPILER_ID STREQUAL "XL")
        set(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
          ${${CMAKE_CXX_COMPILER_ID}_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}})
        string(REPLACE ";" " "
          CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
          "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
        message(STATUS "Setting CXX compiler flags to "
          "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
      else()
        message(STATUS "Unknown CXX compiler ${CMAKE_CXX_COMPILER_ID}. Using default CXX compiler flags")
        set(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
          ${DEFAULT_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}})
        string(REPLACE ";" " "
          CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
          "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
        message(STATUS "Setting CXX compiler flags to "
          "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
      endif()
    else()
      message(STATUS "Using user supplied CXX compiler flags: ${CMAKE_CXX_FLAGS}")
      set(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
        ${CMAKE_CXX_FLAGS})
    endif()
  endif()
endif()

# Additional linker flags
set(MGMOL_LINK_FLAGS "" CACHE STRING "Additional linker flags")
list(APPEND LINK_LIBRARIES ${MGMOL_LINK_FLAGS})

# Use MPI (required)
find_package(MPI REQUIRED)
if (${MPI_CXX_FOUND})
  include_directories(${MPI_INCLUDE_PATH})
  add_definitions(-DUSE_MPI)
else(${MPI_CXX_FOUND})
  message(FATAL_ERROR "Required MPI package not found. Please specify correct path to MPI.")  
endif (${MPI_CXX_FOUND})
message("MPIEXEC :" ${MPIEXEC})
# Use openMP if available
find_package(OpenMP)
if(${OPENMP_CXX_FOUND})
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}
        "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}} -qopenmp")
  endif()
else(${OPENMP_CXX_FOUND})
  message(STATUS "Building without OpenMP")
endif(${OPENMP_CXX_FOUND})

# HDF5 (required)
set(HDF5_PREFER_PARALLEL True)
find_package(HDF5 REQUIRED COMPONENTS C HL)
if (${HDF5_FOUND})
  message(STATUS "HDF5_LIBRARIES after search: ${HDF5_LIBRARIES}")
  message(STATUS "HDF5_HL_LIBRARIES after search: ${HDF5_HL_LIBRARIES}")
else(${HDF5_FOUND})
  message(FATAL_ERROR "Required HDF5 package not found.")
endif (${HDF5_FOUND})

# Boost (required)
find_package(Boost COMPONENTS program_options REQUIRED)
if(NOT ${Boost_FOUND})
  message(FATAL_ERROR "Boost library requited")
endif(NOT ${Boost_FOUND})

# blas/ lapack (required)
# Use cmake option -DBLA_VENDOR=<vendor>
# to set/ find vendor specific libraries.
# (see findBLAS documentation for more info.)
if(NOT BLAS_LIBRARIES)
  find_package(BLAS REQUIRED)
  message(STATUS "BLAS_LIBARIES: ${BLAS_LIBRARIES}")
  if(NOT ${BLAS_FOUND})
    message(FATAL_ERROR "Required blas library not found.")
  endif(NOT ${BLAS_FOUND})  
else(NOT BLAS_LIBRARIES)
  message(STATUS "BLAS_LIBARIES: ${BLAS_LIBRARIES}")
endif(NOT BLAS_LIBRARIES)  

if(NOT LAPACK_LIBRARIES)  
  find_package(LAPACK REQUIRED)
  message(STATUS "LAPACK_LIBARIES: ${LAPACK_LIBRARIES}")
  if(NOT ${LAPACK_FOUND})
    message(FATAL_ERROR "Required LAPACK library not found.")
  endif(NOT ${LAPACK_FOUND}) 
else(NOT LAPACK_LIBRARIES)
  message(STATUS "LAPACK_LIBARIES: ${LAPACK_LIBRARIES}")
endif(NOT LAPACK_LIBRARIES)

# Scalapack (required)
set (MGMOL_WITH_SCALAPACK True CACHE BOOL "Build with SCALAPACK")
if (${MGMOL_WITH_SCALAPACK} OR DEFINED SCALAPACK_ROOT)
  find_package(SCALAPACK)
  if(${SCALAPACK_FOUND})
    add_definitions(-DSCALAPACK)
    message(STATUS "SCALAPACK_INCLUDE_DIRS: ${SCALAPACK_INCLUDE_DIRS}")
    message(STATUS "SCALAPACK_LIBARIES: ${SCALAPACK_LIBRARIES}")
  else(${SCALAPACK_FOUND})
    message(FATAL_ERROR "Required SCALAPACK package not found.")
  endif(${SCALAPACK_FOUND})  
endif(${MGMOL_WITH_SCALAPACK} OR DEFINED SCALAPACK_ROOT)

# ARPACK (optional)
set(MGMOL_WITH_ARPACK FALSE CACHE BOOL "Compile with ARPACK package")
if(${MGMOL_WITH_ARPACK} OR DEFINED ARPACK_ROOT)
  find_package(ARPACK)
  if(${ARPACK_FOUND})
    message(STATUS "Build with ARPACK package")
    add_definitions(-DHAVE_ARPACK)
  endif(${ARPACK_FOUND})
endif(${MGMOL_WITH_ARPACK} OR DEFINED ARPACK_ROOT)

# tricubic (optional)
set(MGMOL_WITH_TRICUBIC FALSE CACHE BOOL "Compile with tricubic package")
if(${MGMOL_WITH_TRICUBIC} OR DEFINED TRICUBIC_ROOT)
  find_package(TRICUBIC)
  if(${TRICUBIC_FOUND})
    message(STATUS "Build with tricubic package")
    add_definitions(-DHAVE_TRICUBIC)
  endif(${TRICUBIC_FOUND})
endif(${MGMOL_WITH_TRICUBIC} OR DEFINED TRICUBIC_ROOT)

enable_testing()

# additional definitions
add_definitions(-DUSE_DIS_MAT)

include(FortranCInterface)
FortranCInterface_HEADER(
    fc_mangle.h
    MACRO_NAMESPACE "FC_"
    SYMBOLS
        DSCAL SSCAL DAXPY SAXPY DDOT SDOT DCOPY SCOPY
        DNRM2 SNRM2 IDAMAX ISAMAX DROT SROT
        DGEMV SGEMV DSYMV
        #blas3
        DGEMM SGEMM DSYMM DSYRK SSYRK DTRMM DTRSM STRSM
        #lapack
        DSYEV DPOTRF DPOTRS DPOTRI DSYGST DTRTRS DPOCON DSYGV DLANGE
)

FortranCInterface_HEADER(
    scalapack_mangle.h
    MACRO_NAMESPACE "FC_"
    SYMBOLS
        numroc indxl2g
        pdswap psswap pdlange pslange pdlacp3 pslacp3 pdtran pstran
        pdsymm pssymm pdgemm pdgemv pdsymv psgemm
        psgemv pssymv pdsyrk pssyrk pdgemr2d psgemr2d pigemr2d
        pdtrmm pstrmm pdtrsm pstrsm pdtrtrs pstrtrs pdpotrf pspotrf
        pdpotrs pspotrs pdgetrf psgetrf pdgetrs psgetrs pdpotri pspotri
        pdtrtri pstrtri pdpocon pspocon pdsygst pssygst pdsyev pssyev
        pdelset pselset pdlatra pslatra pdlaset pslaset pdgesvd psgesvd
        pdamax psamax
)

#FortranCInterface_HEADER(
#    arpack_mangle.h
#    MACRO_NAMESPACE "FC_"
#    SYMBOLS
#        PDSAUPD PDSEUPD DSAUPD DSEUPD
#)

include_directories("${CMAKE_CURRENT_BINARY_DIR}")

execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions(-DGITHASH=${GIT_COMMIT_HASH})

# add include directories
include_directories("${PROJECT_SOURCE_DIR}/src/DistMatrix")
include_directories("${PROJECT_SOURCE_DIR}/src/linear_algebra")
include_directories("${PROJECT_SOURCE_DIR}/src/local_matrices")
include_directories("${PROJECT_SOURCE_DIR}/src/numerical_kernels")
include_directories("${PROJECT_SOURCE_DIR}/src/pb")
include_directories("${PROJECT_SOURCE_DIR}/src/radial")
include_directories("${PROJECT_SOURCE_DIR}/src/sparse_linear_algebra")
include_directories("${PROJECT_SOURCE_DIR}/src/tools")
include_directories("${PROJECT_SOURCE_DIR}/src")

# add subdirectories for source files, tests
add_subdirectory(src)

add_subdirectory(tests)

