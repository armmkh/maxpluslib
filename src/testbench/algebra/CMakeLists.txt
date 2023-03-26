add_test(
    NAME algebra
    COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/testing_algebra
)

add_executable(testing_algebra
    matrixtest.cc
    sparsematrixtest.cc
    testing.cc
    valuetest.cc
    vectortest.cc
)

target_link_libraries(testing_algebra maxplus)

set(MAXPLUSLIB_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MAXPLUSLIB_INCLUDE_DIR}/maxplus
)
