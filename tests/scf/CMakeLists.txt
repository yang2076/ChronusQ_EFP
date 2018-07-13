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

# Directories for SCF tests
set(SCF_TEST_SOURCE_ROOT "${TEST_ROOT}/scf" )
set(SCF_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/scf" )
  

# Set up compilation of SCF test exe
add_executable(scftest ../ut.cxx 
	rhf.cxx uhf.cxx x2chf.cxx ks.cxx 
	rks.cxx uks.cxx x2cks.cxx misc.cxx
        rhfnr.cxx
	rhf_giao.cxx uhf_giao.cxx ghf_giao.cxx)

target_compile_definitions(scftest PUBLIC BOOST_TEST_MODULE=SCF)
target_include_directories(scftest PUBLIC ${SCF_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(scftest PUBLIC ${CQEX_LINK})

if(CQEX_DEP)
  add_dependencies(scftest ${CQEX_DEP})
endif()

# Generate directories
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/rhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/rks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/uhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/uks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/x2c)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/rhf_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/uhf_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/ghf_giao)

file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/rhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/rks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/uhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/uks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/x2c)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/rhf_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/uhf_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/ghf_giao)

# Add the Test
add_test( RHF_SCF   scftest --report_level=detailed --run_test=RHF   )
add_test( UHF_SCF   scftest --report_level=detailed --run_test=UHF   )
add_test( X2CHF_SCF scftest --report_level=detailed --run_test=X2CHF )
add_test( RKS_SCF   scftest --report_level=detailed --run_test=RKS   )
add_test( UKS_SCF   scftest --report_level=detailed --run_test=UKS   )
add_test( X2CKS_SCF scftest --report_level=detailed --run_test=X2CKS )
add_test( MISC_SCF scftest --report_level=detailed --run_test=MISC_SCF)

add_test( RHF_NR_SCF   scftest --report_level=detailed --run_test=RHF_NR   )
       
add_test( RHF_GIAO_SCF   scftest --report_level=detailed --run_test=RHF_GIAO   )
add_test( UHF_GIAO_SCF   scftest --report_level=detailed --run_test=UHF_GIAO   )
add_test( GHF_GIAO_SCF   scftest --report_level=detailed --run_test=GHF_GIAO   )

add_test( KS_KEYWORD scftest --report_level=detailed --run_test=KS_KEYWORD )
add_test( KS_FUNC scftest --report_level=detailed --run_test=KS_FUNC )


if( CQ_ENABLE_MPI )
  add_test( NAME RHF_SCF_MPI   COMMAND  ${MPIEXEC} -np 2 "./scftest" "--report_level=detailed" "--run_test=RHF"     )
  add_test( NAME UHF_SCF_MPI   COMMAND  ${MPIEXEC} -np 2 "./scftest" "--report_level=detailed" "--run_test=UHF"     )
  add_test( NAME X2CHF_SCF_MPI COMMAND  ${MPIEXEC} -np 2 "./scftest" "--report_level=detailed" "--run_test=X2CHF"   )
  #  Seems to be a bug in KS + MPI, need to revisit
  add_test( NAME RKS_SCF_MPI   COMMAND  ${MPIEXEC} -np 2 "./scftest" "--report_level=detailed" "--run_test=RKS"     )
  add_test( NAME UKS_SCF_MPI   COMMAND  ${MPIEXEC} -np 2 "./scftest" "--report_level=detailed" "--run_test=UKS"     )
  add_test( NAME X2CKS_SCF_MPI COMMAND  ${MPIEXEC} -np 2 "./scftest" "--report_level=detailed" "--run_test=X2CKS"   )
endif()
