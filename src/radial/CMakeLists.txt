
set (SOURCES RadialMeshFunction.cc  RadialInter.cc)
add_library(mgmol_radial ${SOURCES})
target_link_libraries(mgmol_radial PUBLIC MPI::MPI_CXX)
install(TARGETS mgmol_radial DESTINATION lib)
