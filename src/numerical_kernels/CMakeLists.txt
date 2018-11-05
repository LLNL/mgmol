
set(SOURCES rho.cc)
add_library(mgmol_numerical_kernels ${SOURCES})
install(TARGETS mgmol_numerical_kernels DESTINATION lib)
