add_test(
    NAME game
    COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/testing_game
)

add_executable(testing_game
    policyiterationtest.cc
    strategyvectortest.cc
    testing.cc
)

target_link_libraries(testing_game maxplus)

set(MAXPLUSLIB_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MAXPLUSLIB_INCLUDE_DIR}/maxplus
)
