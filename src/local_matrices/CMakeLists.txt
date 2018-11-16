
set(SOURCES SquareLocalMatrices.cc LocalMatrices.cc)
add_library(mgmol_local_matrices ${SOURCES}) 
target_compile_definitions (mgmol_local_matrices PRIVATE -DUSE_MPI)
install(TARGETS mgmol_local_matrices DESTINATION lib)
