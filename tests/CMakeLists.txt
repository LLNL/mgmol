find_package(PythonInterp 3.2 REQUIRED)

add_executable(testHDF5P
               ${CMAKE_SOURCE_DIR}/tests/testHDF5P.cc)
add_executable(testMPI
               ${CMAKE_SOURCE_DIR}/tests/testMPI.cc)
add_executable(testBlacsContext
               ${CMAKE_SOURCE_DIR}/tests/testBlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc)
add_executable(testDistVector
               ${CMAKE_SOURCE_DIR}/tests/testDistVector.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc)
add_executable(testDistMatrix
               ${CMAKE_SOURCE_DIR}/tests/testDistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc)
add_executable(testConditionDistMatrix
               ${CMAKE_SOURCE_DIR}/tests/testConditionDistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc)
add_executable(testConditionDistMatrixPower
               ${CMAKE_SOURCE_DIR}/tests/testConditionDistMatrixPower.cc
               ${CMAKE_SOURCE_DIR}/src/Power.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/LocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/SquareLocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/src/tools/random.cc)
add_executable(testPower
               ${CMAKE_SOURCE_DIR}/tests/testPower.cc
               ${CMAKE_SOURCE_DIR}/src/Power.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/LocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/SquareLocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/tools/random.cc)
add_executable(testPowerDistMatrix
               ${CMAKE_SOURCE_DIR}/tests/testPowerDistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/Power.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/LocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/SquareLocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/tools/random.cc)
add_executable(testDirectionalReduce
               ${CMAKE_SOURCE_DIR}/tests/testDirectionalReduce.cc
               ${CMAKE_SOURCE_DIR}/src/sparse_linear_algebra/DirectionalReduce.cc
               ${CMAKE_SOURCE_DIR}/src/pb/PEenv.cc)
add_executable(testAndersonMix
               ${CMAKE_SOURCE_DIR}/tests/Anderson/testAndersonMix.cc
               ${CMAKE_SOURCE_DIR}/tests/Anderson/Solution.cc
               ${CMAKE_SOURCE_DIR}/src/AndersonMix.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc)

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
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testDistVector)
add_test(NAME testDistMatrix
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testDistMatrix)
add_test(NAME testConditionDistMatrix
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testConditionDistMatrix)
add_test(NAME testConditionDistMatrixPower
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testConditionDistMatrixPower)
add_test(NAME testPower
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testPower)
add_test(NAME testPowerDistMatrix
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testPowerDistMatrix)
add_test(NAME testDirectionalReduce
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testDirectionalReduce)
add_test(NAME testAndersonMix
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testAndersonMix 20 2)

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
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/Cl2_ONCVPSP_LDA/mgmol.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/Cl2_ONCVPSP_LDA/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME BandGapN2
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/BandGapN2/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/BandGapN2/mgmol.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/BandGapN2/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/BandGapN2/lrs.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME ReplicatedSP2
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ReplicatedSP2/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/ReplicatedSP2/mgmol.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/ReplicatedSP2/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/ReplicatedSP2/lrs.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testMD_D72
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MD_D72/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 5
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/MD_D72/mgmol_quench.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/MD_D72/mgmol_md.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/MD_D72/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/MD_D72/lrs.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)

target_include_directories(testHDF5P PRIVATE ${HDF5_INCLUDE_DIRS})
target_include_directories(testAndersonMix PRIVATE
                           ${CMAKE_SOURCE_DIR}/tests/Anderson
                           ${HDF5_INCLUDE_DIRS})

target_link_libraries(testHDF5P ${HDF5_LIBRARIES})
target_link_libraries(testBlacsContext ${SCALAPACK_LIBRARIES}
                                       ${BLAS_LIBRARIES})
target_link_libraries(testDistVector ${SCALAPACK_LIBRARIES}
                                     ${BLAS_LIBRARIES})
target_link_libraries(testDistMatrix ${SCALAPACK_LIBRARIES}
                                     ${BLAS_LIBRARIES})
target_link_libraries(testConditionDistMatrix ${SCALAPACK_LIBRARIES}
                                              ${BLAS_LIBRARIES})
target_link_libraries(testConditionDistMatrixPower ${SCALAPACK_LIBRARIES}
                                              ${BLAS_LIBRARIES})
target_link_libraries(testPower ${BLAS_LIBRARIES}
                                ${SCALAPACK_LIBRARIES})
target_link_libraries(testPowerDistMatrix ${BLAS_LIBRARIES}
                                ${SCALAPACK_LIBRARIES})
target_link_libraries(testAndersonMix ${LAPACK_LIBRARIES}
                                      ${BLAS_LIBRARIES}
                                      ${Boost_LIBRARIES})

set_tests_properties(testSiH4 PROPERTIES REQUIRED_FILES
                     ${CMAKE_SOURCE_DIR}/potentials/pseudo.Si)
