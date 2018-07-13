#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2018 Li Research Group (University of Washington)
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

# Directories for basic functionality tests
set(FUNC_TEST_SOURCE_ROOT "${TEST_ROOT}/func" )
set(FUNC_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/func" )


# Set up compilation of Functionality test exe
add_executable(functest ../ut.cxx contract.cxx ordqz.cxx gplhr.cxx)

target_compile_definitions(functest PUBLIC BOOST_TEST_MODULE=FUNC)
target_include_directories(functest PUBLIC ${FUNC_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(functest PUBLIC ${CQEX_LINK})

if(CQEX_DEP)
  add_dependencies(functest ${CQEX_DEP})
endif()

# Add the Tests
add_test( DIRECT_CONTRACTION functest --report_level=detailed --run_test=DIRECT_CONTRACTION)
add_test( ORDQZ functest --report_level=detailed --run_test=ORDQZ)
add_test( GPLHR functest --report_level=detailed --run_test=GPLHR_TESTS)


if( CQ_ENABLE_MPI )

  add_test( NAME GPLHR_MPI COMMAND ${MPIEXEC} -np 2 "./functest" "--report_level=detailed" "--run_test=GPLHR_TESTS" )

endif()
