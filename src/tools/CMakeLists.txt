set(SOURCES mgmol_mpi_tools.cc 
       fermi.cc 
       Vector3D.cc 
       Timer.cc 
       MPIdata.cc 
       MGmol_MPI.cc 
       entropy.cc 
       random.cc 
       coloring.cc 
       SymmetricMatrix.cc
)

add_library(mgmol_tools STATIC ${SOURCES})
target_link_libraries(mgmol_tools PUBLIC MPI::MPI_CXX)
install(TARGETS mgmol_tools DESTINATION lib)
