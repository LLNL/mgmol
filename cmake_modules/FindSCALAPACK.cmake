#---------------
# license file
#---------------
# - FindTricubic module

set(_SCALAPACK_SEARCHES)

# Search SCALAPACK_DIR first if it is set.
if(SCALAPACK_DIR)
  set(_SCALAPACK_SEARCH_DIR PATHS ${SCALAPACK_DIR} NO_DEFAULT_PATH)
  list(APPEND _SCALAPACK_SEARCHES _SCALAPACK_SEARCH_DIR)
endif()

set(SCALAPACK_NAMES scalapack )

# Try each search configuration.
foreach(search ${_SCALAPACK_SEARCHES})
  find_path(SCALAPACK_INCLUDE_DIR NAMES scalapack.h ${${search}} PATH_SUFFIXES include)
endforeach()

# Allow SCALAPACK_LIBRARY to be set manually, as the location of the scalapack library
if(NOT SCALAPACK_LIBRARY)
  foreach(search ${_SCALAPACK_SEARCHES})
    find_library(SCALAPACK_LIBRARY NAMES ${SCALAPACK_NAMES} ${${search}} PATH_SUFFIXES lib)
  endforeach()
endif()

unset(SCALAPACK_NAMES)

mark_as_advanced(SCALAPACK_LIBRARY SCALAPACK_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SCALAPACK REQUIRED_VARS SCALAPACK_LIBRARY)

if(SCALAPACK_FOUND)
  set(SCALAPACK_INCLUDE_DIRS ${SCALAPACK_INCLUDE_DIR})
  
  if(NOT SCALAPACK_LIBRARIES)
    set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY})
  endif(NOT SCALAPACK_LIBRARIES)
endif(SCALAPACK_FOUND)

FIND_LIBRARY(SCALAPACK_BLACS_LIBRARIES 
  NAMES blacsF77init blacsCinit blacs
  PATHS /usr/lib /usr/local/lib /usr/local/tools/lib /usr/local/tools/blacs/lib
  DOC "If a separate BLACS library is required for SCALAPACK, specify it here.")
MARK_AS_ADVANCED(SCALAPACK_BLACS_LIBRARIES)

IF(SCALAPACK_BLACS_LIBRARIES)
  SET(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARIES} ${SCALAPACK_BLACS_LIBRARIES})
ENDIF(SCALAPACK_BLACS_LIBRARIES)
