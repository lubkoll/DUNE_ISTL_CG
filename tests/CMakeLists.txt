project(dune-cg-unit-tests)

cmake_minimum_required(VERSION 2.8)


# Set path to cmake scripts
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    )

# Do not allow in-source builds
include(ConfigSafeGuards)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -coverage -Wall -Wwrite-strings -Wpedantic -std=c++11")
 
# set up include-directories
include_directories(/home/dev/lib/Kaskade_software/gcc-5.2.0/dune-2.4.0/include)
include_directories(${PROJECT_SOURCE_DIR}/../../../FGlue-c++11)
include_directories(${PROJECT_SOURCE_DIR}/..)

# google test
set( GTEST_DIR /home/dev/Dokumente/Paper/2016/Lubkoll-FunG/experiments/FunG/external )
include_directories(${GTEST_DIR}/googletest ${GTEST_DIR}/googletest/include )

add_custom_target(git_update
                COMMAND git submodule init
                COMMAND git submodule update
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
add_library(gtest
            ${GTEST_DIR}/googletest/src/gtest-all.cc
            ${GTEST_DIR}/googletest/src/gtest_main.cc)

add_dependencies(gtest git_update)

set_source_files_properties(${GTEST_DIR}/googletest/src/gtest-all.cc  PROPERTIES GENERATED 1)
set_source_files_properties(${GTEST_DIR}/googletest/src/gtest_main.cc PROPERTIES GENERATED 1)

aux_source_directory(. SRC_LIST)
aux_source_directory(mock SRC_LIST)

add_executable(unit_tests ${SRC_LIST})
target_link_libraries(unit_tests gtest pthread)

include(CTest)
enable_testing()
add_test(test ${PROJECT_BINARY_DIR}/unit_tests)