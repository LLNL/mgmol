
add_subdirectory(DistMatrix)
add_subdirectory(linear_algebra)
add_subdirectory(local_matrices)
add_subdirectory(numerical_kernels)
add_subdirectory(pb)
add_subdirectory(radial)
add_subdirectory(sparse_linear_algebra)
add_subdirectory(tools)

set(link_libs mgmol_distmatrix 
   mgmol_linear_algebra 
   mgmol_local_matrices 
   mgmol_numerical_kernels 
   mgmol_pb 
   mgmol_radial 
   mgmol_sparse_linear_algebra 
   mgmol_tools
   )

set(SOURCES
 LocalMatrices2DistMatrix.cc
 DistMatrix2SquareLocalMatrices.cc
 GrassmanCGFactory.cc
 DMStrategyFactory.cc
 manage_memory.cc
 SubCell.cc
 SinCosOps.cc 
 PolakRibiereSolver.cc 
 EnergySpreadPenalty.cc 
 SpreadPenalty.cc 
 SpreadPenaltyVolume.cc 
 setup.cc 
 MVP_DMStrategy.cc 
 HamiltonianMVP_DMStrategy.cc 
 MVPSolver.cc 
 HamiltonianMVPSolver.cc 
 OrbitalsPreconditioning.cc 
 DFTsolver.cc 
 NonOrthoDMStrategy.cc 
 FullyOccupiedNonOrthoDMStrategy.cc 
 EigenDMStrategy.cc 
 Masks4Orbitals.cc 
 AOMMprojector.cc 
 hdf_tools.cc 
 DistributedIonicData.cc 
 IonData.cc 
 SubspaceProjector.cc 
 MDfiles.cc 
 OrbitalsExtrapolation.cc 
 OrbitalsExtrapolationOrder2.cc 
 OrbitalsExtrapolationOrder3.cc 
 runfire.cc 
 FIRE.cc 
 IonicAlgorithm.cc 
 FIRE_IonicStepper.cc 
 tools.cc 
 MGmol.cc 
 MGmol_NEB.cc 
 ABPG.cc 
 GrassmanLineMinimization.cc 
 GrassmanCG.cc 
 GrassmanCGSparse.cc 
 LBFGS.cc 
 IonicStepper.cc 
 Energy.cc 
 GramMatrix.cc 
 DensityMatrix.cc 
 DensityMatrixSparse.cc 
 ConstraintSet.cc 
 DistanceConstraint.cc 
 MultiDistanceConstraint.cc 
 getLAeigen.cc 
 LDAFunctional.cc 
 LDAonGrid.cc 
 PBEFunctional.cc 
 PBEonGrid.cc 
 PBEonGridSpin.cc 
 Electrostatic.cc 
 ProjectedMatrices.cc 
 ProjectedMatricesSparse.cc 
 ProjectedMatrices2N.cc 
 ProjectedMatricesMehrstellen.cc 
 ReplicatedWorkSpace.cc 
 Rho.cc 
 Mesh.cc 
 Control.cc 
 Hamiltonian.cc 
 BlockVector.cc 
 KBPsiMatrixInterface.cc 
 KBPsiMatrixSparse.cc 
 KBprojectorSparse.cc 
 MultipoleExpansion.cc 
 SpreadsAndCenters.cc 
 Preconditioning.cc 
 OrbitalsTransform.cc 
 NOLMOTransform.cc 
 LocalizationRegions.cc 
 Hartree.cc 
 ShiftedHartree.cc 
 PBdiel.cc 
 Species.cc 
 HDFrestart.cc 
 MasksSet.cc 
 FunctionsPacking.cc 
 ColoredRegions.cc 
 LocGridOrbitals.cc 
 ExtendedGridOrbitals.cc
 mcstep.cc 
 MD_IonicStepper.cc 
 LBFGS_IonicStepper.cc 
 lbfgsrlx.cc 
 OrthoAndersonMix.cc
 AndersonMix.cc 
 MLWFTransform.cc 
 Ion.cc 
 GridMask.cc 
 GridMaskMult.cc 
 GridMaskMax.cc 
 Ions.cc 
 restart.cc 
 md.cc 
 get_vnlpsi.cc 
 quench.cc 
 mlwf.cc 
 readInput.cc 
 Forces.cc 
 computeHij.cc 
 Potentials.cc 
 BasicDataDistributors.cc 
 ShortSightedInverse.cc 
 PCGSolver_Diel.cc 
 PBdiel_CG.cc 
 PCGSolver.cc 
 Hartree_CG.cc 
 ClusterOrbitals.cc
 SP2.cc
 Power.cc
 PowerGen.cc
 SuperSampling.cc
 magma_singleton.cc
)

add_library(mgmol_src ${SOURCES})

target_include_directories(mgmol_src PRIVATE ${HDF5_INCLUDE_DIRS})
target_include_directories(mgmol_src PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories (mgmol_src PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(mgmol_src ${link_libs})
if(${MGMOL_WITH_MAGMA})
  target_link_libraries(mgmol_src PkgConfig::MAGMA)
endif()

# compile definitions
#target_compile_definitions (mgmol_src PUBLIC -DUSE_CNR_)
target_compile_definitions (mgmol_src PUBLIC -DUSE_MPI)

install(TARGETS mgmol_src DESTINATION lib)

# build executable
add_executable(mgmol-opt main.cc)
target_include_directories (mgmol-opt PRIVATE ${Boost_INCLUDE_DIRS})

target_link_libraries(mgmol-opt mgmol_src ${link_libs}) 
target_link_libraries(mgmol-opt ${SCALAPACK_LIBRARIES})
target_link_libraries(mgmol-opt ${HDF5_LIBRARIES})
target_link_libraries(mgmol-opt ${HDF5_HL_LIBRARIES})
target_link_libraries(mgmol-opt ${BLAS_LIBRARIES})
target_link_libraries(mgmol-opt ${LAPACK_LIBRARIES})
target_link_libraries(mgmol-opt ${Boost_LIBRARIES})
if (${OPENMP_CXX_FOUND})
  target_link_libraries(mgmol-opt OpenMP::OpenMP_CXX)
endif()

install(TARGETS mgmol-opt DESTINATION bin)

