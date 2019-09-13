# - Find the MAGMA library
#
# Usage:
#   find_package(MAGMA [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   MAGMA_FOUND               ... true if magma is found on the system
#   MAGMA_LIBRARY_DIRS        ... full path to magma library
#   MAGMA_INCLUDE_DIRS        ... magma include directory
#   MAGMA_LIBRARIES           ... magma libraries


find_package(PkgConfig REQUIRED)
pkg_check_modules(MAGMA REQUIRED magma IMPORTED_TARGET)
