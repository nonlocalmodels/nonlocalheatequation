#  Copyright (c) 2014 Hartmut Kaiser
#  Copyright (c) 2014 Patricia Grubel
#  Copyright (c) 2020 Pranav Gadikar
#
#  SPDX-License-Identifier: BSL-1.0
#  Distributed under the Boost Software License, Version 1.0. (See accompanying
#  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


PROJECT(Nonlocalheatequation)

cmake_minimum_required(VERSION 3.10)

SET(BUILD_PROJECT_DESCRIPTION OFF CACHE STRING "Build project description")

if(BUILD_PROJECT_DESCRIPTION)
    add_subdirectory(description)
endif()

find_package(HPX REQUIRED)
add_hpx_executable(1d_nonlocal_serial
    SOURCES ./src/1d_nonlocal_serial.cpp
    COMPONENT_DEPENDENCIES iostreams)