
set(SOURCES mputils.cc  serial_functions.cc)
add_library(mgmol_linear_algebra STATIC ${SOURCES})
target_link_libraries(mgmol_linear_algebra PUBLIC MPI::MPI_CXX)
install(TARGETS mgmol_linear_algebra DESTINATION lib)
