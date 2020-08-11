#  Copyright (c) 2014 Hartmut Kaiser
#  Copyright (c) 2014 Patricia Grubel
#  Copyright (c) 2020 Pranav Gadikar
#
#  SPDX-License-Identifier: BSL-1.0
#  Distributed under the Boost Software License, Version 1.0. (See accompanying
#  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


PROJECT(Nonlocalheatequation)
set (VERSION_MAJOR 0)
set (VERSION_MINOR 1)
set (VERSION_UPDATE 0)

configure_file (
  "${PROJECT_SOURCE_DIR}/Config.h.in"
  "${PROJECT_SOURCE_DIR}/include/Config.h"
  )

cmake_minimum_required(VERSION 3.10)

set(BUILD_PROJECT_DESCRIPTION OFF CACHE STRING "Build project description")
set(CMAKE_BUILD_TYPE RelWithDebInfo)

if(BUILD_PROJECT_DESCRIPTION)
    add_subdirectory(description)
endif()

find_package(HPX REQUIRED)
find_package(VTK REQUIRED)

# Set directory for all exectuables
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")

# Set the common library directory
get_property(LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS)

if ("${LIB64}" STREQUAL "TRUE")
    set(LIBSUFFIX 64)
else()
    set(LIBSUFFIX "")
endif()

set(LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib${LIBSUFFIX}")

# Add a postfix if the libraries are build in the debug mode
set(CMAKE_DEBUG_POSTFIX d)

# user needs to set GMSH_DIR and METIS_DIR variables if the libraries are installed
# in custom locations
find_library(GMSH_LIBRARY NAMES gmsh
    PATHS /usr
    /usr/local
    ${GMSH_DIR}
    PATH_SUFFIXES lib lib64)
find_library(METIS_LIBRARY NAMES metis
    PATHS /usr
    /usr/local
    ${METIS_DIR}
    PATH_SUFFIXES lib lib64)

include(${VTK_USE_FILE})

add_library(LibsModule 
    ./include/writer.cpp
)
target_link_libraries(LibsModule ${VTK_LIBRARIES})

add_executable(2d_domain_decomposition ./src/domain_decomposition.cpp)
target_link_libraries(2d_domain_decomposition ${METIS_LIBRARY} ${GMSH_LIBRARY})

add_hpx_executable(1d_nonlocal_serial
    SOURCES ./src/1d_nonlocal_serial.cpp
    DEPENDENCIES LibsModule
    COMPONENT_DEPENDENCIES iostreams)

add_hpx_executable(2d_nonlocal_serial
    SOURCES ./src/2d_nonlocal_serial.cpp
    DEPENDENCIES LibsModule
    COMPONENT_DEPENDENCIES iostreams)

add_hpx_executable(2d_nonlocal_async
    SOURCES ./src/2d_nonlocal_async.cpp
    DEPENDENCIES LibsModule
    COMPONENT_DEPENDENCIES iostreams)

add_hpx_executable(2d_nonlocal_distributed
    SOURCES ./src/2d_nonlocal_distributed.cpp
    DEPENDENCIES LibsModule
    COMPONENT_DEPENDENCIES iostreams)

enable_testing()

# test for checking the correctness of 1d nonlocal equation's serial implementation
# test what         : correctness for large epsilon, small epsilon, larger number timesteps, various time and space discretizations
# test how          : parameters will be inputted from the file '1d.txt' file in 'tests' folder as described below
# first line        : number of test cases
# from second line  : actual parameters inputs to be tested
# order of inputs   : nx, nt, eps, k, dt, dx
# Note              : Please refer to 'program_options' in 'main' function for '1d_nonlocal_serial.cpp' file in 'src' folder
add_test(NAME Test_1d COMMAND sh -c "./bin/1d_nonlocal_serial --test_batch < ../tests/1d.txt")
set_tests_properties(Test_1d
    PROPERTIES
    PASS_REGULAR_EXPRESSION "Tests Passed")
set_tests_properties(Test_1d
    PROPERTIES
    FAIL_REGULAR_EXPRESSION "Tests Failed")

# test for checking the correctness of 2d nonlocal equation's serial implementation
# test what         : correctness for large epsilon, small epsilon, larger number timesteps, various time and space discretizations
# test how          : parameters will be inputted from the file '2d.txt' file in 'tests' folder as described below
# first line        : number of test cases
# from second line  : actual parameters inputs to be tested
# order of inputs   : nx, ny, nt, eps, k, dt, dh
# Note              : Please refer to 'program_options' in 'main' function for '2d_nonlocal_serial.cpp' file in 'src' folder
add_test(NAME Test_2d COMMAND sh -c "./bin/2d_nonlocal_serial --test_batch < ../tests/2d.txt")
set_tests_properties(Test_2d
    PROPERTIES
    PASS_REGULAR_EXPRESSION "Tests Passed")
set_tests_properties(Test_2d
    PROPERTIES
    FAIL_REGULAR_EXPRESSION "Tests Failed")

# test for checking the correctness of 2d nonlocal equation's asynchronous implementation
# test what         : correctness for large epsilon, small epsilon, large number paritions, small number of paritions
#                     larger number timesteps, various time and space discretizations
# test how          : parameters will be inputted from the file '2d_async.txt' file in 'tests' folder as described below
# first line        : number of test cases
# from second line  : actual parameters inputs to be tested
# order of inputs   : nx, ny, np, nt, eps, k, dt, dh
# Note              : Please refer to 'program_options' in 'main' function for '2d_nonlocal_async.cpp' file in 'src' folder
add_test(NAME Test_2d_async COMMAND sh -c "./bin/2d_nonlocal_async --test_batch < ../tests/2d_async.txt")
set_tests_properties(Test_2d_async
    PROPERTIES
    PASS_REGULAR_EXPRESSION "Tests Passed")
set_tests_properties(Test_2d_async
    PROPERTIES
    FAIL_REGULAR_EXPRESSION "Tests Failed")

# test for checking the correctness of 2d distributed equation's asynchronous implementation
# test what         : correctness for large epsilon, small epsilon, large number paritions, small number of paritions
#                     larger number timesteps, various time and space discretizations
# test how          : parameters will be inputted from the file '2d_distributed.txt' file in 'tests' folder as described below
# first line        : number of test cases
# from second line  : actual parameters inputs to be tested
# order of inputs   : nx, ny, npx, npy, nt, eps, k, dt, dh
# Note              : Please refer to 'program_options' in 'main' function for '2d_nonlocal_async.cpp' file in 'src' folder
add_test(NAME Test_2d_distributed COMMAND sh -c "./bin/2d_nonlocal_distributed --test_batch < ../tests/2d_distributed.txt")
set_tests_properties(Test_2d_distributed
    PROPERTIES
    PASS_REGULAR_EXPRESSION "Tests Passed")
set_tests_properties(Test_2d_distributed
    PROPERTIES
    FAIL_REGULAR_EXPRESSION "Tests Failed")

# Install the executables and libs
install(TARGETS 2d_nonlocal_distributed 
        LIBRARY DESTINATION lib${LIBSUFFIX}
        ARCHIVE DESTINATION lib${LIBSUFFIX}
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)
