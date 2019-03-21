
set(SOURCES SquareLocalMatrices.cc LocalMatrices.cc)
add_library(mgmol_local_matrices STATIC ${SOURCES}) 
install(TARGETS mgmol_local_matrices DESTINATION lib)
