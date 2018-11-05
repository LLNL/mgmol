#---------------
# license file
#---------------
# - FindARPACK module

set(_ARPACK_SEARCHES)

# Search ARPACK_DIR first if it is set.
if(ARPACK_DIR)
  set(_ARPACK_SEARCH_DIR PATHS ${ARPACK_DIR} NO_DEFAULT_PATH)
  list(APPEND _ARPACK_SEARCHES _ARPACK_SEARCH_DIR)
endif()

set(ARPACK_NAMES parpack arpack )

# Try each search configuration.
foreach(search ${_ARPACK_SEARCHES})
  find_path(ARPACK_INCLUDE_DIR NAMES ARPACK.h ${${search}} PATH_SUFFIXES include)
endforeach()

# Allow ARPACK_LIBRARY to be set manually, as the location of the ARPACK library
if(NOT ARPACK_LIBRARY)
  foreach(search ${_ARPACK_SEARCHES})
    find_library(ARPACK_LIBRARY NAMES ${ARPACK_NAMES} ${${search}} PATH_SUFFIXES lib)
  endforeach()
endif()

unset(ARPACK_NAMES)

mark_as_advanced(ARPACK_LIBRARY ARPACK_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ARPACK REQUIRED_VARS ARPACK_LIBRARY)

if(ARPACK_FOUND)
  set(ARPACK_INCLUDE_DIRS ${ARPACK_INCLUDE_DIR})
  
  if(NOT ARPACK_LIBRARIES)
    set(ARPACK_LIBRARIES ${ARPACK_LIBRARY})
  endif(NOT ARPACK_LIBRARIES)
endif(ARPACK_FOUND)
