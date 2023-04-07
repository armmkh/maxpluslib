add_test(
    NAME base
    COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/testing_base
)

add_executable(testing_base
    testing.cc
    mcmtest.cc
)

target_link_libraries(testing_base maxplus)

set(MAXPLUSLIB_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MAXPLUSLIB_INCLUDE_DIR}/maxplus
)
