set(_SCALAPACK_SEARCHES)

# Search SCALAPACK_ROOT first if it is set.
if(DEFINED ENV{SCALAPACK_ROOT})
  set(_SCALAPACK_SEARCH_DIR $ENV{SCALAPACK_ROOT} $ENV{SCALAPACK_ROOT}/lib/intel64)
  list(APPEND _SCALAPACK_SEARCHES ${_SCALAPACK_SEARCH_DIR})
endif(DEFINED ENV{SCALAPACK_ROOT})

if(SCALAPACK_ROOT)
  set(_SCALAPACK_SEARCH_DIR ${SCALAPACK_ROOT} ${SCALAPACK_ROOT}/lib/intel64)
  list(APPEND _SCALAPACK_SEARCHES ${_SCALAPACK_SEARCH_DIR})
endif(SCALAPACK_ROOT)

set(SCALAPACK_NAMES scalapack mkl_scalapack_lp64 scalapack-openmpi scalapack-pvm
  scalapack-mpi scalapack-mpich scalapack-mpich2 scalapack-lam)

# Try each search configuration.
if(NOT SCALAPACK_INCLUDE_DIR)
  foreach(search ${_SCALAPACK_SEARCHES})
    find_path(SCALAPACK_INCLUDE_DIR NAMES scalapack.h mkl_scalapack.h PATHS ${search} PATH_SUFFIXES include NO_DEFAULT_PATH)
    message(STATUS "SCALAPACK_INCLUDE_DIR: ${SCALAPACK_INCLUDE_DIR}")
  endforeach()
endif(NOT SCALAPACK_INCLUDE_DIR)

# Allow SCALAPACK_LIBRARY to be set manually, as the location of the scalapack library
if(NOT SCALAPACK_LIBRARY)
  foreach(search ${_SCALAPACK_SEARCHES})
    find_library(SCALAPACK_LIBRARY NAMES ${SCALAPACK_NAMES} PATHS ${search} PATH_SUFFIXES lib NO_DEFAULT_PATH)
  endforeach()
endif()

mark_as_advanced(SCALAPACK_LIBRARY SCALAPACK_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SCALAPACK REQUIRED_VARS SCALAPACK_LIBRARY)

# Search for some default library paths
if (NOT SCALAPACK_FOUND)
  find_library(SCALAPACK_LIBRARY
    NAMES ${SCALAPACK_NAMES}
    PATHS /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib
    /opt/local/lib /opt/sw/lib /sw/lib
    ENV LD_LIBRARY_PATH
    ENV DYLD_FALLBACK_LIBRARY_PATH
    ENV DYLD_LIBRARY_PATH
    ENV SCALAPACKDIR
    ENV BLACSDIR)

  FIND_PACKAGE_HANDLE_STANDARD_ARGS(SCALAPACK REQUIRED_VARS SCALAPACK_LIBRARY)
endif()

unset(SCALAPACK_NAMES)

if(SCALAPACK_FOUND)
  # Only Intel's scalapack requires an include directory
  if(SCALAPACK_INCLUDE_DIR)
    set(SCALAPACK_INCLUDE_DIRS ${SCALAPACK_INCLUDE_DIR})
  endif(SCALAPACK_INCLUDE_DIR)
  
  if(NOT SCALAPACK_LIBRARIES)
    set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY})
  endif(NOT SCALAPACK_LIBRARIES)
endif(SCALAPACK_FOUND)

FIND_LIBRARY(SCALAPACK_BLACS_LIBRARY 
  NAMES blacsF77init blacsCinit blacs
  PATHS /usr/lib /usr/local/lib /usr/local/tools/lib /usr/local/tools/blacs/lib
  DOC "If a separate BLACS library is required for SCALAPACK, specify it here.")
MARK_AS_ADVANCED(SCALAPACK_BLACS_LIBRARY)

IF(SCALAPACK_BLACS_LIBRARY)
  SET(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARIES} ${SCALAPACK_BLACS_LIBRARY})
ENDIF(SCALAPACK_BLACS_LIBRARY)
