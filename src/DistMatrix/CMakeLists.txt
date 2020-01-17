
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
target_include_directories (mgmol_distmatrix PRIVATE ${SCALAPACK_INCLUDE_DIRS})
target_link_libraries(mgmol_distmatrix PUBLIC MPI::MPI_CXX)
if(${MGMOL_WITH_MAGMA})
  target_link_libraries(mgmol_distmatrix PRIVATE PkgConfig::MAGMA)
endif()

target_compile_definitions (mgmol_distmatrix PUBLIC -DSCALAPACK)
install(TARGETS mgmol_distmatrix DESTINATION lib) 
