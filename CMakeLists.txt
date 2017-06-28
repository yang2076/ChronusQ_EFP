#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2017 Li Research Group (University of Washington)
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 
# Contact the Developers:
#   E-Mail: xsli@uw.edu
#
cmake_minimum_required(VERSION 3.0)
project(ChronusQ CXX C)
set(ChronusQ_VERSION_MAJOR 0)
set(ChronusQ_VERSION_MINOR 1)
set(ChronusQ_VERSION_PATCH 0)


# Append local cmake directory to find CMAKE Modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")


# C++11
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Define include directories
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")

include(FindLibint)
include(Boost)
include(HandleLA)
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CQEX_LINK ${CQ_EXT_LINK} )



set(BASIS_PATH     "${PROJECT_SOURCE_DIR}/basis"     )
set(AUTHOR_PATH    "${PROJECT_SOURCE_DIR}/AUTHORS"   )

configure_file(
  "${PROJECT_SOURCE_DIR}/include/chronusq_config.hpp.in"
  "${PROJECT_BINARY_DIR}/include/chronusq_config.hpp"
)

# Recurse into subdirectories
add_subdirectory(src)


