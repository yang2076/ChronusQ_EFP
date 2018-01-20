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
include(ExternalProject)


find_package(Eigen3 REQUIRED)
message(STATUS "Found Eigen3 -- ${EIGEN3_INCLUDE_DIR}")
include_directories("${EIGEN3_INCLUDE_DIR}")

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  message(STATUS "INTEL COMPILER -- Adding MKL to Compiler Invocation")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=parallel -std=c++11")
  set(_CQ_MKL 1)
else()

  set( OPENBLAS_PREFIX ${PROJECT_SOURCE_DIR}/external/openblas )
  set( OPENBLAS_INCLUDEDIR ${OPENBLAS_PREFIX}/include )
  set( OPENBLAS_LIBDIR ${OPENBLAS_PREFIX}/lib )

  if( NOT EXISTS "${OPENBLAS_INCLUDEDIR}/f77blas.h" )


    if( OPENBLAS_TARGET )
      message( STATUS "Forcing OpenBLAS TARGET = ${OPENBLAS_TARGET}" )
      set(OPENBLAS_BUILD_COMMAND make TARGET=${OPENBLAS_TARGET} -j2)
    else()
      message( STATUS "Allowing OpenBLAS to determine CPU TARGET" )
      set(OPENBLAS_BUILD_COMMAND make -j2)
    endif()

    ExternalProject_Add(openblas
      PREFIX ${OPENBLAS_PREFIX}
      URL "${OPENBLAS_PREFIX}/v0.2.19.tar.gz"
      CONFIGURE_COMMAND echo 'No OpenBLAS Configure Command'
      BUILD_COMMAND ${OPENBLAS_BUILD_COMMAND}
      BUILD_IN_SOURCE 1
      INSTALL_COMMAND make install PREFIX=${OPENBLAS_PREFIX} &&
        cd ${OPENBLAS_INCLUDEDIR} && patch < ../patch/lapacke.patch &&
        patch < ../patch/f77blas.patch 
    )


    list(APPEND CQEX_DEP openblas)
  endif()

  include_directories(${OPENBLAS_INCLUDEDIR})
  link_directories(${OPENBLAS_LIBDIR})
  list(APPEND CQ_EXT_LINK ${OPENBLAS_LIBDIR}/libopenblas.a gfortran)
endif()
