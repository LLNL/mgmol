
set(SOURCES rho.cc)
add_library(mgmol_numerical_kernels ${SOURCES})
target_link_libraries(mgmol_numerical_kernels PUBLIC MPI::MPI_CXX)
install(TARGETS mgmol_numerical_kernels DESTINATION lib)
