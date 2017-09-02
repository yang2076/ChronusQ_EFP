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

# Directories for RT tests
set(RT_TEST_SOURCE_ROOT "${TEST_ROOT}/rt" )
set(RT_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/rt" )
  

# Set up compilation of RT test exe
add_executable(rttest ../ut.cxx rrt.cxx urt.cxx grt.cxx)

target_compile_definitions(rttest PUBLIC BOOST_TEST_MODULE=RT)
target_include_directories(rttest PUBLIC ${RT_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(rttest PUBLIC ${CQEX_LINK})

if(CQEX_DEP)
  add_dependencies(rttest ${CQEX_DEP})
endif()

# Generate directories
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/serial/rrt)
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/serial/urt)
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/serial/grt)

file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/parallel/rrt)
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/parallel/urt)
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/parallel/grt)

# Add the Test
add_test( RHF_RT   rttest --report_level=detailed --run_test=RHF_RT   )
add_test( UHF_RT   rttest --report_level=detailed --run_test=UHF_RT   )
add_test( RKS_RT   rttest --report_level=detailed --run_test=RKS_RT   )
add_test( UKS_RT   rttest --report_level=detailed --run_test=UKS_RT   )
add_test( X2CHF_RT rttest --report_level=detailed --run_test=X2CHF_RT )
