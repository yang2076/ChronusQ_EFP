include(ExternalProject)


find_package(Eigen3)
include_directories("${EIGEN3_INCLUDE_DIR}")

if(${CMAKE_CXX_COMPILER} MATCHES "icpc.*$")
  message(STATUS "INTEL COMPILER")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=parallel -std=c++11")
  set(_CQ_MKL 1)
else()

  set( OPENBLAS_PREFIX ${PROJECT_SOURCE_DIR}/external/openblas )
  set( OPENBLAS_INCLUDEDIR ${OPENBLAS_PREFIX}/include )
  set( OPENBLAS_LIBDIR ${OPENBLAS_PREFIX}/lib )

  if( NOT EXISTS "${OPENBLAS_INCLUDEDIR}/f77blas.h" )

    ExternalProject_Add(openblas
      PREFIX ${OPENBLAS_PREFIX}
      URL "${OPENBLAS_PREFIX}/v0.2.19.tar.gz"
      CONFIGURE_COMMAND echo 'No OpenBLAS Configure Command'
      BUILD_COMMAND make -j2
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
