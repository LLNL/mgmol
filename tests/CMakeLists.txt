find_package(PythonInterp 3.2 REQUIRED)

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
               ${CMAKE_SOURCE_DIR}/src/tools/random.cc
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
               ${CMAKE_SOURCE_DIR}/src/tools/random.cc
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
               ${CMAKE_SOURCE_DIR}/src/tools/random.cc
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
               ${CMAKE_SOURCE_DIR}/src/tools/random.cc
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
               ${CMAKE_SOURCE_DIR}/tests/ut_serial_main.cc)
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
               ${CMAKE_SOURCE_DIR}/tests/ut_serial_main.cc)
add_executable(testVariableSizeMatrix
               ${CMAKE_SOURCE_DIR}/tests/testVariableSizeMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/sparse_linear_algebra/VariableSizeMatrix.cc               
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/sparse_linear_algebra/SparseRow.cc               
               ${CMAKE_SOURCE_DIR}/src/sparse_linear_algebra/SparseRowAndTable.cc     
               ${CMAKE_SOURCE_DIR}/src/sparse_linear_algebra/Table.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testTradeGhostValues
               ${CMAKE_SOURCE_DIR}/tests/testTradeGhostValues.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/pb/Grid.cc
               ${CMAKE_SOURCE_DIR}/src/pb/PEenv.cc
               ${CMAKE_SOURCE_DIR}/src/pb/GridFunc.cc
               ${CMAKE_SOURCE_DIR}/src/pb/GridFuncVector.cc
               ${CMAKE_SOURCE_DIR}/src/pb/MGkernels.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testLaph4
               ${CMAKE_SOURCE_DIR}/tests/testLaph4.cc
               ${CMAKE_SOURCE_DIR}/src/pb/Laph4.cc
               ${CMAKE_SOURCE_DIR}/src/pb/Lap.cc
               ${CMAKE_SOURCE_DIR}/src/pb/Laph2.cc
               ${CMAKE_SOURCE_DIR}/src/pb/FDoper.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/pb/Grid.cc
               ${CMAKE_SOURCE_DIR}/src/pb/PEenv.cc
               ${CMAKE_SOURCE_DIR}/src/pb/GridFunc.cc
               ${CMAKE_SOURCE_DIR}/src/pb/GridFuncVector.cc
               ${CMAKE_SOURCE_DIR}/src/pb/MGkernels.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testBatchLaph4
               ${CMAKE_SOURCE_DIR}/tests/testBatchLaph4.cc
               ${CMAKE_SOURCE_DIR}/src/pb/Laph4.cc
               ${CMAKE_SOURCE_DIR}/src/pb/Lap.cc
               ${CMAKE_SOURCE_DIR}/src/pb/Laph2.cc
               ${CMAKE_SOURCE_DIR}/src/pb/FDoper.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/pb/Grid.cc
               ${CMAKE_SOURCE_DIR}/src/pb/PEenv.cc
               ${CMAKE_SOURCE_DIR}/src/pb/GridFunc.cc
               ${CMAKE_SOURCE_DIR}/src/pb/GridFuncVector.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/src/pb/MGkernels.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testMGkernels
               ${CMAKE_SOURCE_DIR}/tests/testMGkernels.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/pb/Grid.cc
               ${CMAKE_SOURCE_DIR}/src/pb/PEenv.cc
               ${CMAKE_SOURCE_DIR}/src/pb/GridFunc.cc
               ${CMAKE_SOURCE_DIR}/src/pb/GridFuncVector.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/src/pb/MGkernels.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
add_executable(testGramMatrix
               ${CMAKE_SOURCE_DIR}/tests/testGramMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/GramMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/Power.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/LocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/SquareLocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/ReplicatedMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/MatricesBlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrixTools.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/src/tools/random.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_magma_main.cc)
add_executable(testDensityMatrix
               ${CMAKE_SOURCE_DIR}/tests/testDensityMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DensityMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/GramMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/Power.cc
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/LocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/local_matrices/SquareLocalMatrices.cc
               ${CMAKE_SOURCE_DIR}/src/ReplicatedMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrix.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/BlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/MatricesBlacsContext.cc
               ${CMAKE_SOURCE_DIR}/src/DistMatrix/DistMatrixTools.cc
               ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
               ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
               ${CMAKE_SOURCE_DIR}/src/tools/MGmol_MPI.cc
               ${CMAKE_SOURCE_DIR}/src/tools/mgmol_mpi_tools.cc
               ${CMAKE_SOURCE_DIR}/src/tools/random.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_magma_main.cc)

if(${MAGMA_FOUND})
  add_executable(testOpenmpOffload
               ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
               ${CMAKE_SOURCE_DIR}/tests/testOpenmpOffload.cc
               ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
  if(${MGMOL_WITH_OPENMP_OFFLOAD})
    add_executable(testMagmaOpenmp
                   ${CMAKE_SOURCE_DIR}/tests/testMagmaOpenmp.cc
                   ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
  endif()
  add_executable(testReplicatedMatrix
                 ${CMAKE_SOURCE_DIR}/tests/testReplicatedMatrix.cc
                 ${CMAKE_SOURCE_DIR}/src/ReplicatedMatrix.cc
                 ${CMAKE_SOURCE_DIR}/src/local_matrices/LocalMatrices.cc
                 ${CMAKE_SOURCE_DIR}/src/local_matrices/SquareLocalMatrices.cc
                 ${CMAKE_SOURCE_DIR}/src/magma_singleton.cc
                 ${CMAKE_SOURCE_DIR}/src/linear_algebra/mputils.cc
                 ${CMAKE_SOURCE_DIR}/src/tools/Timer.cc
                 ${CMAKE_SOURCE_DIR}/src/tools/random.cc
                 ${CMAKE_SOURCE_DIR}/tests/ut_main.cc)
endif()
               
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
add_test(NAME testVariableSizeMatrix
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testVariableSizeMatrix)
add_test(NAME testTradeGhostValues
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testTradeGhostValues)
add_test(NAME testLaph4
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testLaph4)
add_test(NAME testBatchLaph4
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testBatchLaph4)
add_test(NAME testtMGkernels
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testMGkernels)
add_test(NAME testGramMatrix
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testGramMatrix)
add_test(NAME testDensityMatrix
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testDensityMatrix)

if(${MAGMA_FOUND}) 
  add_test(NAME testOpenmpOffload
           COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testOpenmpOffload)
  if(${MGMOL_WITH_OPENMP_OFFLOAD})
    add_test(NAME testMagmaOpenmp
             COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testMagmaOpenmp)
  endif()
  add_test(NAME testReplicatedMatrix
           COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                   ${CMAKE_CURRENT_BINARY_DIR}/testReplicatedMatrix)
else()
  add_test(NAME testShortSighted
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ShortSighted/test.py
           ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 5 ${MPIEXEC_PREFLAGS}
           ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
           ${CMAKE_CURRENT_SOURCE_DIR}/ShortSighted/quench.cfg
           ${CMAKE_CURRENT_SOURCE_DIR}/ShortSighted/md.cfg
           ${CMAKE_CURRENT_SOURCE_DIR}/ShortSighted/coords.in
           ${CMAKE_CURRENT_SOURCE_DIR}/ShortSighted/lrs.in
           ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
  add_test(NAME testProjectedMatrices
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ProjectedMatrices/test.py
           ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
           ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
           ${CMAKE_CURRENT_SOURCE_DIR}/ProjectedMatrices/quenchExact.cfg
           ${CMAKE_CURRENT_SOURCE_DIR}/ProjectedMatrices/quenchShortSighted.cfg
           ${CMAKE_CURRENT_SOURCE_DIR}/ProjectedMatrices/coords.in
           ${CMAKE_CURRENT_SOURCE_DIR}/ProjectedMatrices/lrs.in
           ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
endif()

#Regression tests
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
add_test(NAME testSpreadPenalty
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpreadPenalty/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/SpreadPenalty/mgmol.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/SpreadPenalty/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/SpreadPenalty/lrs.in
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
add_test(NAME testLBFGS
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/LBFGS/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/LBFGS/lbfgs.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/LBFGS/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/LBFGS/lrs.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME ChebyshevMVP
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Chebyshev/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/Chebyshev/cheb.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/Chebyshev/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testContinuumSolvent
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ContinuumSolvent/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/ContinuumSolvent/mgmol_diel.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/ContinuumSolvent/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/ContinuumSolvent/lrs.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testSpinO2
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpinO2/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/SpinO2/mgmol_spin1.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/SpinO2/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)
add_test(NAME testSpinO2LDA
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SpinO2LDA/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
         ${CMAKE_CURRENT_BINARY_DIR}/../src/mgmol-opt
         ${CMAKE_CURRENT_SOURCE_DIR}/SpinO2LDA/mgmol.cfg
         ${CMAKE_CURRENT_SOURCE_DIR}/SpinO2LDA/coords.in
         ${CMAKE_CURRENT_SOURCE_DIR}/SpinO2LDA/lrs.in
         ${CMAKE_CURRENT_SOURCE_DIR}/../potentials)


target_include_directories(testHDF5P PRIVATE ${HDF5_INCLUDE_DIRS})
target_include_directories(testAndersonMix PRIVATE
                           ${CMAKE_SOURCE_DIR}/tests/Anderson
                           ${HDF5_INCLUDE_DIRS})

target_link_libraries(testHDF5P PRIVATE ${HDF5_LIBRARIES}
  MPI::MPI_CXX)
target_link_libraries(testMPI PRIVATE MPI::MPI_CXX)
target_link_libraries(testBlacsContext PRIVATE ${SCALAPACK_LIBRARIES}
  ${BLAS_LIBRARIES} MPI::MPI_CXX)
target_link_libraries(testSuperSampling PRIVATE MPI::MPI_CXX)
target_link_libraries(testDirectionalReduce PRIVATE MPI::MPI_CXX)

if(${MAGMA_FOUND})
  target_link_libraries(testDistVector PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testReplicated2DistMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testDistMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testConditionDistMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testConditionDistMatrixPower PRIVATE
    ${SCALAPACK_LIBRARIES} ${BLAS_LIBRARIES} 
    MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testReplicatedMatrix PRIVATE
    MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testPower PRIVATE ${BLAS_LIBRARIES} ${SCALAPACK_LIBRARIES}
    MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testPowerDistMatrix PRIVATE ${BLAS_LIBRARIES}
    ${SCALAPACK_LIBRARIES} MPI::MPI_CXX
    OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testAndersonMix PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} 
    MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testVariableSizeMatrix PRIVATE ${BLAS_LIBRARIES} 
    MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testTradeGhostValues PRIVATE ${BLAS_LIBRARIES}
    MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testLaph4 PRIVATE ${BLAS_LIBRARIES}
    MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testBatchLaph4 PRIVATE ${BLAS_LIBRARIES}
    MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testMGkernels PRIVATE ${BLAS_LIBRARIES}
    MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testGramMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testDensityMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  target_link_libraries(testOpenmpOffload PRIVATE 
    MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)

  if(${MGMOL_WITH_OPENMP_OFFLOAD})
    target_link_libraries(testMagmaOpenmp PRIVATE 
      MPI::MPI_CXX OpenMP::OpenMP_CXX PkgConfig::MAGMA)
  endif()
else()
  target_link_libraries(testDistVector PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testReplicated2DistMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testDistMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testConditionDistMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testConditionDistMatrixPower PRIVATE
    ${SCALAPACK_LIBRARIES} ${BLAS_LIBRARIES} 
    MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testPower PRIVATE ${BLAS_LIBRARIES} ${SCALAPACK_LIBRARIES}
    MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testPowerDistMatrix PRIVATE ${BLAS_LIBRARIES}
    ${SCALAPACK_LIBRARIES} MPI::MPI_CXX
    OpenMP::OpenMP_CXX)
  target_link_libraries(testAndersonMix PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} 
    MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testVariableSizeMatrix PRIVATE ${BLAS_LIBRARIES} 
    MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testTradeGhostValues PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
    MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testLaph4 PRIVATE ${BLAS_LIBRARIES}
    MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testBatchLaph4 PRIVATE ${BLAS_LIBRARIES}
    MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testMGkernels PRIVATE ${BLAS_LIBRARIES}
    MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testGramMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX)
  target_link_libraries(testDensityMatrix PRIVATE ${SCALAPACK_LIBRARIES}
    ${BLAS_LIBRARIES} MPI::MPI_CXX OpenMP::OpenMP_CXX)
endif()

set_tests_properties(testSiH4 PROPERTIES REQUIRED_FILES
                     ${CMAKE_SOURCE_DIR}/potentials/pseudo.Si)
