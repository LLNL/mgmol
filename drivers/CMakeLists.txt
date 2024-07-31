include_directories( ${CMAKE_SOURCE_DIR}/src )

add_executable(check_input check_input.cc)

add_executable(example1 example1.cc)

target_include_directories(check_input PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(example1 PRIVATE ${Boost_INCLUDE_DIRS})
if(DEFINED LBFGSPP_INCLUDE)
  add_executable(lbfgspp lbfgspp.cc)

  target_include_directories(lbfgspp PRIVATE ${Boost_INCLUDE_DIRS})
  target_include_directories(lbfgspp PRIVATE ${EIGEN_INCLUDE})
  target_include_directories(lbfgspp PRIVATE ${LBFGSPP_INCLUDE})
endif(DEFINED LBFGSPP_INCLUDE)

target_link_libraries(check_input mgmol_src)
target_link_libraries(example1 mgmol_src)
target_link_libraries(lbfgspp mgmol_src)
