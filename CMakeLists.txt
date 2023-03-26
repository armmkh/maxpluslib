cmake_minimum_required(VERSION 3.14.0 FATAL_ERROR)
project(MaxPLusLib)

option(CODE_COVERAGE "Compile for code coverage (default OFF)." OFF)

set(CPM_USE_LOCAL_PACKAGES ON)
include(config/get_cpm.cmake)

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
  endif()

  # Let's ensure -std=c++xx instead of -std=g++xx
  set(CMAKE_CXX_EXTENSIONS OFF)
  set(CMAKE_CXX_STANDARD 17)

  # # Set output folders
  set(MAXPLUSLIB_BUILD_DIRECTORY ${CMAKE_BINARY_DIR}/build) 
  set(MAXPLUSLIB_COVERAGE_DIRECTORY ${MAXPLUSLIB_BUILD_DIRECTORY}/coverage) 
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${MAXPLUSLIB_BUILD_DIRECTORY}/lib) 
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${MAXPLUSLIB_BUILD_DIRECTORY}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${MAXPLUSLIB_BUILD_DIRECTORY}/bin)

  # Let's nicely support folders in IDEs
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)

  # Testing only available if this is the main app
  # Note this needs to be done in the main CMakeLists
  # since it calls enable_testing, which must be in the
  # main CMakeLists.
  include(CTest)

  find_package(Doxygen OPTIONAL_COMPONENTS dot)
  if (Doxygen_FOUND)
    add_subdirectory(documentation)
  else()
    message(STATUS "Doxygen not found, not building docs")
  endif()


endif()

### CODE COVERAGE ###
if(${CODE_COVERAGE})
  message("Compiling for code coverage analysis.")
  set(CMAKE_BUILD_TYPE Debug)

  include(config/CodeCoverage.cmake)
  append_coverage_compiler_flags()
  # add_compile_options("-fprofile-arcs")
  # add_compile_options("-ftest-coverage")
endif()

enable_testing()

add_custom_target(maxpluslibcoverage_setup
  ${CMAKE_COMMAND} -E make_directory ${MAXPLUSLIB_COVERAGE_DIRECTORY}
)

add_custom_target(maxpluslibcoverage1
  lcov  --capture --directory ${CMAKE_BINARY_DIR} --output ${MAXPLUSLIB_BUILD_DIRECTORY}/coverage.info
  DEPENDS maxpluslibcoverage_setup
)

add_custom_target(maxpluslibcoverage2
  lcov -r ${MAXPLUSLIB_BUILD_DIRECTORY}/coverage.info /usr/include/\* -o ${MAXPLUSLIB_BUILD_DIRECTORY}/coverage.info
  DEPENDS maxpluslibcoverage1
)

add_custom_target(maxpluslibcoverage
  genhtml ${MAXPLUSLIB_BUILD_DIRECTORY}/coverage.info --output-directory ${MAXPLUSLIB_COVERAGE_DIRECTORY}/
  DEPENDS maxpluslibcoverage2
)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(src)
