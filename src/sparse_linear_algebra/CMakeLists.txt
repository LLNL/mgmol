
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

target_include_directories(mgmol_sparse_linear_algebra PUBLIC ${HDF5_INCLUDE_DIRS})

install(TARGETS mgmol_sparse_linear_algebra DESTINATION lib)
