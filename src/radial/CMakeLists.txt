
set (SOURCES RadialFunction.cc  RadialInter.cc)
add_library(mgmol_radial ${SOURCES})
install(TARGETS mgmol_radial DESTINATION lib)
