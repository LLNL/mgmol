
set (SOURCES RadialFunction.cc  RadialInter.cc)
add_library(mgmol_radial STATIC ${SOURCES})
target_compile_definitions (mgmol_radial PRIVATE -DUSE_MPI)
install(TARGETS mgmol_radial DESTINATION lib)
