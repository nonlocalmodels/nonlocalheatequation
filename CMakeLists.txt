PROJECT(Nonlocalheatequation)

cmake_minimum_required(VERSION 3.16)

SET(BUILD_PROJECT_DESCRIPTION OFF CACHE STRING "Build project description")


if(BUILD_PROJECT_DESCRIPTION)
    add_subdirectory(description)
endif()
