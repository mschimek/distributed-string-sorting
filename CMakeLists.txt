################################################################################
# CMakeLists.txt
#
# Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/misc/cmake")

project(dsss)
add_compile_options(-std=c++17)
add_compile_options(-O3)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(DSSS_OVERSUBSCRIBE
  "Run tests with mpirun [...] --oversubscribe.
  This may be necessary for environments with less than four cores." OFF)

# needs to be set on fedora?
set(CMAKE_PREFIX_PATH /usr/lib64/openmpi)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

# TODO: Look for a "modern" (target based) working way
find_package(MPI REQUIRED)
if (MPI_FOUND)
  message(STATUS "MPI FOUND: ${MPI_C_INCLUDE_PATH}")
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  message("blabla " ${MPI_LIBRARIES})
endif()
#add_subdirectory(external/dsss)
#add_subdirectory(external/benchmark)
#add_subdirectory(external/external)
#add_subdirectory(external/tests)
include_directories(src)
include_directories(tlx)
add_subdirectory(tlx)
add_subdirectory(src)

message(STATUS "Built Type: " ${CMAKE_BUILD_TYPE} )

################################################################################
