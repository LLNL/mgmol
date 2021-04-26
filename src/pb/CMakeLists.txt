
set(SOURCES GridFuncVector.cc 
       Grid.cc 
       GridFunc.cc 
       TriCubic.cc 
       PB.cc 
       SolverPB.cc 
       SolverLap.cc 
       tools.cc 
       PEenv.cc 
       DielFunc.cc 
       FDoper.cc 
       Lap.cc 
       Laph2.cc 
       Laph4.cc 
       Laph4M.cc 
       Laph6.cc 
       Laph8.cc 
       ShiftedLaph4M.cc 
       PBh2.cc 
       PBh4.cc 
       PBh4M.cc 
       PBh6.cc 
       PBh8.cc
       MGkernels.cc
       FDkernels.cc
)

add_library(mgmol_pb ${SOURCES})
target_include_directories (mgmol_pb PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mgmol_pb PUBLIC MPI::MPI_CXX)
target_link_libraries(mgmol_pb PUBLIC OpenMP::OpenMP_CXX)
if(${MGMOL_WITH_MAGMA})
  target_link_libraries(mgmol_pb PRIVATE PkgConfig::MAGMA)
endif()
install(TARGETS mgmol_pb DESTINATION lib)
