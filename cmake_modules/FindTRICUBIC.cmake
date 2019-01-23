#---------------
# license file
#---------------
# - FindTricubic module

set(_TRICUBIC_SEARCHES)

# Search TRICUBIC_ROOT first if it is set.
if(TRICUBIC_ROOT)
  set(_TRICUBIC_SEARCH_DIR PATHS ${TRICUBIC_ROOT} NO_DEFAULT_PATH)
  list(APPEND _TRICUBIC_SEARCHES _TRICUBIC_SEARCH_DIR)
endif()

set(TRICUBIC_NAMES tricubic )

# Try each search configuration.
foreach(search ${_TRICUBIC_SEARCHES})
  find_path(TRICUBIC_INCLUDE_DIR NAMES tricubic.h ${${search}} PATH_SUFFIXES include)
endforeach()

# Allow TRICUBIC_LIBRARY to be set manually, as the location of the tricubic library
if(NOT TRICUBIC_LIBRARY)
  foreach(search ${_TRICUBIC_SEARCHES})
    find_library(TRICUBIC_LIBRARY NAMES ${TRICUBIC_NAMES} ${${search}} PATH_SUFFIXES lib)
  endforeach()
endif()

unset(TRICUBIC_NAMES)

mark_as_advanced(TRICUBIC_LIBRARY TRICUBIC_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(TRICUBIC REQUIRED_VARS TRICUBIC_LIBRARY)

if(TRICUBIC_FOUND)
  set(TRICUBIC_INCLUDE_DIRS ${TRICUBIC_INCLUDE_DIR})
  
  if(NOT TRICUBIC_LIBRARIES)
    set(TRICUBIC_LIBRARIES ${TRICUBIC_LIBRARY})
  endif(NOT TRICUBIC_LIBRARIES)
endif(TRICUBIC_FOUND)
