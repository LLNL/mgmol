
set(SOURCES MPI_AllToSome.cc 
       BlacsContext.cc 
       DistMatrix.cc 
       SubMatricesIndexing.cc 
       SubMatrices.cc 
       SparseDistMatrix.cc 
       MatricesBlacsContext.cc 
       DistMatrixTools.cc
  )
  
add_library(mgmol_distmatrix ${SOURCES})
target_include_directories (mgmol_distmatrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions (mgmol_distmatrix PUBLIC -DUSE_DIS_MAT)
install(TARGETS mgmol_distmatrix DESTINATION lib) 
