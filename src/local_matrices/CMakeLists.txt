
set(SOURCES SquareLocalMatrices.cc LocalMatrices.cc SquareSubMatrix.cc)
add_library(mgmol_local_matrices ${SOURCES}) 
target_link_libraries(mgmol_local_matrices PUBLIC MPI::MPI_CXX)
target_link_libraries(mgmol_local_matrices PUBLIC OpenMP::OpenMP_CXX)
if(${MGMOL_WITH_MAGMA})
  target_link_libraries(mgmol_local_matrices PRIVATE PkgConfig::MAGMA)
endif()
install(TARGETS mgmol_local_matrices DESTINATION lib)
