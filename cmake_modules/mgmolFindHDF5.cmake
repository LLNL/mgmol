#---------------
# license file
#---------------
# - mgmolFindHDF5 module

set(HDF5_FOUND FALSE)

set(_HDF5_SEARCHES)

if(NOT HDF5_LIBRARIES)
  if(DEFINED ENV{HDF5_ROOT})
    # Search HDF5_ROOT first if it is set.
    set(_HDF5_SEARCH_DIR $ENV{HDF5_ROOT})
    list(APPEND _HDF5_SEARCHES ${_HDF5_SEARCH_DIR})

    set(HDF5_NAMES "hdf5;hdf5_hl")

    # Try each search configuration.
    if(NOT HDF5_INCLUDE_DIR)
#      file(GLOB HEADER_FILES *.h)
      foreach(search ${_HDF5_SEARCHES})
        find_path(HDF5_INCLUDE_DIR NAMES hdf5.h PATHS ${search} PATH_SUFFIXES include NO_DEFAULT_PATH)
      endforeach()
    endif(NOT HDF5_INCLUDE_DIR)
    # Allow HDF5_LIBRARY to be set manually, as the location of the hdf5 library
    if(NOT HDF5_LIBRARY)
      set(HDF5_LIBRARY)
      foreach(libname ${HDF5_NAMES})
        foreach(search ${_HDF5_SEARCHES})
          find_library(HDF5_LIB_${libname} NAMES ${libname} PATHS ${search} PATH_SUFFIXES lib NO_DEFAULT_PATH)
        endforeach()
        list(APPEND HDF5_LIBRARY ${HDF5_LIB_${libname}})
      endforeach()
    endif()

    unset(HDF5_NAMES)

    mark_as_advanced(HDF5_LIBRARY HDF5_INCLUDE_DIR)

    include(FindPackageHandleStandardArgs)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(HDF5 REQUIRED_VARS HDF5_LIBRARY HDF5_INCLUDE_DIR)

    if(HDF5_FOUND)
      set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
  
      if(NOT HDF5_LIBRARIES)
        set(HDF5_LIBRARIES ${HDF5_LIBRARY})
        set(HDF5_CXX_LIBRARIES ${HDF5_LIBRARY})
      endif(NOT HDF5_LIBRARIES)
    endif(HDF5_FOUND)  

  else(DEFINED ENV{HDF5_ROOT})
    set(HDF5_PREFER_PARALLEL True)
    set(MGMOL_HDF5_COMPONENTS CXX HL)
    find_package(HDF5 COMPONENTS ${MGMOL_HDF5_COMPONENTS} REQUIRED)
    ## find_package will potentially find multiple libraries.
    ## Ensure that a correct library is found.
    if(HDF5_LIBRARIES AND NOT HDF5_IS_PARALLEL)
      set(HDF5_FOUND FALSE)
    endif(HDF5_LIBRARIES AND NOT HDF5_IS_PARALLEL)

  endif(DEFINED ENV{HDF5_ROOT})
  
  # Additional libraries: add libz
  find_library(LIBZ NAMES z PATH_SUFFIXES lib)
  if(LIBZ)
    SET(HDF5_LIBRARIES ${HDF5_LIBRARIES} ${LIBZ})
    SET(HDF5_CXX_LIBRARIES ${HDF5_CXX_LIBRARIES} ${LIBZ})
  endif(LIBZ)
endif(NOT HDF5_LIBRARIES)

