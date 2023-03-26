add_test(
    NAME graph
    COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/testing_graph
)

add_executable(testing_graph
    matrixlabeledfsmtest.cc
    mpautomatontest.cc
    testing.cc
)

target_link_libraries(testing_graph maxplus)

set(MAXPLUSLIB_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MAXPLUSLIB_INCLUDE_DIR}/maxplus
)