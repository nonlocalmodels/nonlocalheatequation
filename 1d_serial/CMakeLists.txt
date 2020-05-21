cmake_minimum_required(VERSION 3.3.2)
project(sample CXX)
find_package(HPX REQUIRED)
add_hpx_executable(1d_nonlocal_serial
    SOURCES 1d_nonlocal_serial.cpp
    COMPONENT_DEPENDENCIES iostreams)
