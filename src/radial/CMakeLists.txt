
set (SOURCES RadialMeshFunction.cc  RadialInter.cc)
add_library(mgmol_radial STATIC ${SOURCES})
install(TARGETS mgmol_radial DESTINATION lib)
