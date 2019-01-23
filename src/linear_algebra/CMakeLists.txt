
set(SOURCES mputils.cc  serial_functions.cc)
add_library(mgmol_linear_algebra STATIC ${SOURCES})
install(TARGETS mgmol_linear_algebra DESTINATION lib)
