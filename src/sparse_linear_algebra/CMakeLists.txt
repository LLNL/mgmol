
set(SOURCES DirectionalReduce.cc 
       SparseRow.cc 
       SparseRowAndTable.cc 
       Table.cc 
       PackedCommunicationBuffer.cc 
       DataDistribution.cc 
       VariableSizeMatrix.cc 
       LinearSolverMatrix.cc 
       PreconILU.cc 
       LinearSolver.cc
)

add_library(mgmol_sparse_linear_algebra STATIC ${SOURCES})

if(${HDF5_FOUND})
  target_include_directories(mgmol_sparse_linear_algebra PUBLIC ${HDF5_INCLUDE_DIRS})
endif(${HDF5_FOUND})

target_compile_definitions (mgmol_sparse_linear_algebra PRIVATE -DUSE_MPI)
install(TARGETS mgmol_sparse_linear_algebra DESTINATION lib)
