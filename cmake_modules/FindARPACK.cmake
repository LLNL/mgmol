#---------------
# license file
#---------------
# - FindARPACK module

set(_ARPACK_SEARCHES)

# Search ARPACK_ROOT first if it is set.
if(ARPACK_ROOT)
  set(_ARPACK_SEARCH_DIR PATHS ${ARPACK_ROOT} NO_DEFAULT_PATH)
  list(APPEND _ARPACK_SEARCHES _ARPACK_SEARCH_DIR)
endif()

set(ARPACK_NAMES parpack arpack )

# Allow ARPACK_LIBRARY to be set manually, as the location of the ARPACK library
if(NOT ARPACK_LIBRARY)
  foreach(search ${_ARPACK_SEARCHES})
    find_library(ARPACK_LIBRARY NAMES ${ARPACK_NAMES} ${${search}} PATH_SUFFIXES lib)
  endforeach()
endif()

unset(ARPACK_NAMES)

mark_as_advanced(ARPACK_LIBRARY)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ARPACK REQUIRED_VARS ARPACK_LIBRARY)

if(ARPACK_FOUND)
  if(NOT ARPACK_LIBRARIES)
    set(ARPACK_LIBRARIES ${ARPACK_LIBRARY})
  endif(NOT ARPACK_LIBRARIES)
endif(ARPACK_FOUND)
