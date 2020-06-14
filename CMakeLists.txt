#  Copyright (c) 2014 Hartmut Kaiser
#  Copyright (c) 2014 Patricia Grubel
#  Copyright (c) 2020 Pranav Gadikar
#
#  SPDX-License-Identifier: BSL-1.0
#  Distributed under the Boost Software License, Version 1.0. (See accompanying
#  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


PROJECT(Nonlocalheatequation)

cmake_minimum_required(VERSION 3.10)

set(BUILD_PROJECT_DESCRIPTION OFF CACHE STRING "Build project description")
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(VTK_DIR "/home/pranav/vtk/build")

if(BUILD_PROJECT_DESCRIPTION)
    add_subdirectory(description)
endif()

find_package(HPX REQUIRED)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_library(LibsModule 
    ./include/writer.cpp
)
target_link_libraries(LibsModule ${VTK_LIBRARIES})

add_hpx_executable(1d_nonlocal_serial
    SOURCES ./src/1d_nonlocal_serial.cpp
    COMPONENT_DEPENDENCIES iostreams)

add_hpx_executable(2d_nonlocal_serial
    SOURCES ./src/2d_nonlocal_serial.cpp
    COMPONENT_DEPENDENCIES iostreams)

enable_testing()

add_test(NAME Test_1d COMMAND sh -c "./1d_nonlocal_serial --test_batch < ../tests/1d.txt")
set_tests_properties(Test_1d
    PROPERTIES
    PASS_REGULAR_EXPRESSION "Tests Passed")
set_tests_properties(Test_1d
    PROPERTIES
    FAIL_REGULAR_EXPRESSION "Tests Failed")

add_test(NAME Test_2d COMMAND sh -c "./2d_nonlocal_serial --test_batch < ../tests/2d.txt")
set_tests_properties(Test_2d
    PROPERTIES
    PASS_REGULAR_EXPRESSION "Tests Passed")
set_tests_properties(Test_2d
    PROPERTIES
    FAIL_REGULAR_EXPRESSION "Tests Failed")