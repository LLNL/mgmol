
set(SOURCES SquareLocalMatrices.cc LocalMatrices.cc)
add_library(mgmol_local_matrices STATIC ${SOURCES}) 
target_link_libraries(mgmol_local_matrices PUBLIC MPI::MPI_CXX)
if(${MGMOL_WITH_MAGMA})
  target_link_libraries(mgmol_local_matrices PRIVATE PkgConfig::MAGMA)
endif()
install(TARGETS mgmol_local_matrices DESTINATION lib)
