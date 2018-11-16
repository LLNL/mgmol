
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

set(SOURCES  SubCell.cc 
 SinCosOps.cc 
 Power.cc 
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
 KBPsiMatrix.cc 
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
 mcstep.cc 
 MD_IonicStepper.cc 
 LBFGS_IonicStepper.cc 
 lbfgsrlx.cc 
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
 DavidsonSolver.cc
 ChebyshevApproximation.cc
 ChebyshevApproximationSparse.cc
 OrbitalsExtrapolationReversible.cc
 LocalizationRegionsXLBOMD.cc
 AuxiliaryLocGridOrbitals.cc
 SubspaceAlignment.cc
)

add_library(mgmol_src ${SOURCES})

if(${HDF5_FOUND})
  target_include_directories(mgmol_src PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(mgmol_src ${HDF5_LIBRARIES})
endif(${HDF5_FOUND})

if(${Boost_FOUND})
  target_include_directories(mgmol_src PUBLIC ${Boost_INCLUDE_DIRS})
  target_link_libraries(mgmol_src ${Boost_LIBRARIES})
endif(${Boost_FOUND})

if(${BLAS_FOUND})
  target_link_libraries(mgmol_src ${BLAS_LIBRARIES})
  target_link_libraries(mgmol_src ${LAPACK_LIBRARIES})
endif(${BLAS_FOUND})

target_include_directories (mgmol_src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mgmol_src ${link_libs})

# compile definitions
target_compile_definitions (mgmol_src PRIVATE -D__MGMOL_DAVIDSON__)
target_compile_definitions (mgmol_src PRIVATE -D__MGMOL_CHEBYSHEV__)
target_compile_definitions (mgmol_src PRIVATE -D__MGMOL_XLBOMD__)
target_compile_definitions (mgmol_src PRIVATE -DHAVE_BOOST)
target_compile_definitions (mgmol_src PRIVATE -DADD_)
target_compile_definitions (mgmol_src PRIVATE -DUSE_CNR_)
target_compile_definitions (mgmol_src PRIVATE -DUSE_MPI)
# compile options
#target_compile_options (mgmol_src PRIVATE -Wl --no_undefined)

#target_include_directories (mgmol_src PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/DistMatrix)

#target_include_directories(mgmol_src PUBLIC DistMatrix)
#add_dependencies(mgmol_src mgmol_distmatrix)
install(TARGETS mgmol_src DESTINATION lib)

   

add_executable(mgmol-opt main.cc)
target_link_libraries(mgmol-opt mgmol_src ${link_libs}) 

install(TARGETS mgmol-opt DESTINATION bin)

