
set (SOURCES RadialFunction.cc  RadialInter.cc)
add_library(mgmol_radial ${SOURCES})
target_compile_definitions (mgmol_radial PRIVATE -DUSE_MPI)
install(TARGETS mgmol_radial DESTINATION lib)
