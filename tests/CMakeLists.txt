find_package(PythonInterp 3.2 REQUIRED)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

add_executable(testHDF5P
               ${CMAKE_SOURCE_DIR}/tests/testHDF5P.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testMPI
               ${CMAKE_SOURCE_DIR}/tests/testMPI.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testBlacsContext
               ${CMAKE_SOURCE_DIR}/tests/testBlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testDistVector
               ${CMAKE_SOURCE_DIR}/tests/testDistVector.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testReplicated2DistMatrix
               ${CMAKE_SOURCE_DIR}/tests/testReplicated2DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/LocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/SquareLocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testDistMatrix
               ${CMAKE_SOURCE_DIR}/tests/testDistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testConditionDistMatrix
               ${CMAKE_SOURCE_DIR}/tests/testConditionDistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testConditionDistMatrixPower
               ${CMAKE_SOURCE_DIR}/tests/testConditionDistMatrixPower.cc
               ${CMAKE_SOURCE_DIR}/src/Power.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/LocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/SquareLocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/src/tools/random.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testPower
               ${CMAKE_SOURCE_DIR}/tests/testPower.cc
               ${CMAKE_SOURCE_DIR}/src/Power.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/LocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/SquareLocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/tools/random.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testPowerDistMatrix
               ${CMAKE_SOURCE_DIR}/tests/testPowerDistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/Power.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/LocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/SquareLocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/tools/random.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testDirectionalReduce
               ${CMAKE_SOURCE_DIR}/tests/testDirectionalReduce.cc
               ${CMAKE_SOURCE_DIR}/src/sparse_linear_algebra/DirectionalReduce.cc
               ${CMAKE_SOURCE_DIR}/src/pb/PEenv.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testAndersonMix
               ${CMAKE_SOURCE_DIR}/tests/Anderson/testAndersonMix.cc
               ${CMAKE_SOURCE_DIR}/tests/Anderson/Solution.cc
               ${CMAKE_SOURCE_DIR}/src/AndersonMix.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc)
add_executable(testSuperSampling
               ${CMAKE_SOURCE_DIR}/tests/testSuperSampling.cc
               ${CMAKE_SOURCE_DIR}/src/SuperSampling.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testVariableSizeMatrix
               ${CMAKE_SOURCE_DIR}/tests/testVariableSizeMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/sparse_linear_algebra/VariableSizeMatrix.cc               
               ${CMAKE_SOURCE_DIR}/src/sparse_linear_algebra/SparseRow.cc               
               ${CMAKE_SOURCE_DIR}/src/sparse_linear_algebra/SparseRowAndTable.cc     
               ${CMAKE_SOURCE_DIR}/src/sparse_linear_algebra/Table.cc                                        
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
               
target_compile_definitions(testAndersonMix PUBLIC TESTING)

add_test(NAME testHDF5P
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testHDF5P)
add_test(NAME testMPI
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testMPI)
add_test(NAME testBlacsContext
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testBlacsContext)
add_test(NAME testDistVector
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testDistVector)
add_test(NAME testReplicated2DistMatrix
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testReplicated2DistMatrix)
add_test(NAME testDistMatrix
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testDistMatrix)
add_test(NAME testConditionDistMatrix
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testConditionDistMatrix)
add_test(NAME testConditionDistMatrixPower
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testConditionDistMatrixPower)
add_test(NAME testPower
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testPower)
add_test(NAME testPowerDistMatrix
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testPowerDistMatrix)
add_test(NAME testDirectionalReduce
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testDirectionalReduce)
add_test(NAME testAndersonMix
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testAndersonMix 20 2)
add_test(NAME testSuperSampling
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testSuperSampling)

 add_test(NAME testProjectedMatrices
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ProjectedMatrices/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/ProjectedMatrices/quenchExact.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/ProjectedMatrices/quenchShortSighted.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/ProjectedMatrices/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/ProjectedMatrices/lrs.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME Davidson
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Davidson/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/Davidson/davidson.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/Davidson/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testMVP
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MVP/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/MVP/mvp.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/MVP/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testFatom
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Fatom/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/Fatom/mgmol.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/Fatom/F.xyz
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testSiH4
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SiH4/testSiH4.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/SiH4/mgmol.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/SiH4/sih4.xyz
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testCl2
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Cl2_ONCVPSP_LDA/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/Cl2_ONCVPSP_LDA/mgmol.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/Cl2_ONCVPSP_LDA/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME BandGapN2
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/BandGapN2/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/BandGapN2/mgmol.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/BandGapN2/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/BandGapN2/lrs.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME ReplicatedSP2
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ReplicatedSP2/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/ReplicatedSP2/mgmol.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/ReplicatedSP2/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/ReplicatedSP2/lrs.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testMLWF
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MLWF/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/MLWF/mlwf.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/MLWF/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testMD_D72
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MD_D72/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 5 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/MD_D72/mgmol_quench.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/MD_D72/mgmol_md.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/MD_D72/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/MD_D72/lrs.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testShortSighted
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ShortSighted/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 5 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/ShortSighted/quench.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/ShortSighted/md.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/ShortSighted/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/ShortSighted/lrs.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testLBFGS
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/LBFGS/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/LBFGS/lbfgs.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/LBFGS/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/LBFGS/lrs.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testVariableSizeMatrix
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testVariableSizeMatrix)

target_include_directories(testHDF5P PRIVATE ${HDF5_INCLUDE_DIRS})
target_include_directories(testAndersonMix PRIVATE
                           ${CMAKE_SOURCE_DIR}/tests/Anderson
                           ${HDF5_INCLUDE_DIRS})

target_link_libraries(testHDF5P PRIVATE ${HDF5_LIBRARIES}
  Boost::unit_test_framework MPI::MPI_CXX)
target_link_libraries(testMPI PRIVATE Boost::unit_test_framework MPI::MPI_CXX)
target_link_libraries(testBlacsContext PRIVATE ${SCALAPACK_LIBRARIES}
  ${BLAS_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX)
target_link_libraries(testSuperSampling PRIVATE Boost::unit_test_framework MPI::MPI_CXX)
target_link_libraries(testDirectionalReduce PRIVATE Boost::unit_test_framework MPI::MPI_CXX)
target_link_libraries(testVariableSizeMatrix PRIVATE ${BLAS_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX)

if(${MAGMA_FOUND})
  target_link_libraries(testDistVector PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX PkgConfig::MAGMA)
  target_link_libraries(testReplicated2DistMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX PkgConfig::MAGMA)
  target_link_libraries(testDistMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX PkgConfig::MAGMA)
  target_link_libraries(testConditionDistMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX PkgConfig::MAGMA)
  target_link_libraries(testConditionDistMatrixPower PRIVATE
    ${SCALAPACK_LIBRARIES} ${BLAS_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX PkgConfig::MAGMA)
  target_link_libraries(testPower PRIVATE ${BLAS_LIBRARIES} ${SCALAPACK_LIBRARIES}
    Boost::unit_test_framework MPI::MPI_CXX PkgConfig::MAGMA)
  target_link_libraries(testPowerDistMatrix PRIVATE ${BLAS_LIBRARIES}
    ${SCALAPACK_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX
    PkgConfig::MAGMA)
  target_link_libraries(testAndersonMix PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} 
    ${Boost_LIBRARIES} MPI::MPI_CXX PkgConfig::MAGMA)
else()
  target_link_libraries(testDistVector PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX)
  target_link_libraries(testReplicated2DistMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX)
  target_link_libraries(testDistMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX)
  target_link_libraries(testConditionDistMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX)
  target_link_libraries(testConditionDistMatrixPower PRIVATE
    ${SCALAPACK_LIBRARIES} ${BLAS_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX)
  target_link_libraries(testPower PRIVATE ${BLAS_LIBRARIES} ${SCALAPACK_LIBRARIES}
    Boost::unit_test_framework MPI::MPI_CXX)
  target_link_libraries(testPowerDistMatrix PRIVATE ${BLAS_LIBRARIES}
    ${SCALAPACK_LIBRARIES} Boost::unit_test_framework MPI::MPI_CXX)
  target_link_libraries(testAndersonMix PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} 
    ${Boost_LIBRARIES} MPI::MPI_CXX)
endif()

set_tests_properties(testSiH4 PROPERTIES REQUIRED_FILES
                     ${CMAKE_SOURCE_DIR}/potentials/pseudo.Si)
