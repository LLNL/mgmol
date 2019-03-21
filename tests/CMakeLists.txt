
add_executable(testDirectionalReduce
               ${CMAKE_SOURCE_DIR}/tests/testDirectionalReduce.cc
               ${CMAKE_SOURCE_DIR}/src/sparse_linear_algebra/DirectionalReduce.cc
               ${CMAKE_SOURCE_DIR}/src/pb/PEenv.cc)

add_test(NAME testDirectionalReduce
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testDirectionalReduce)

