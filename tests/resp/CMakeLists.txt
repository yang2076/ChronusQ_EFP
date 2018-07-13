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

# Directories for RESP tests
set(RESP_TEST_SOURCE_ROOT "${TEST_ROOT}/resp" )
set(RESP_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/resp" )
  

# Set up compilation of RESP test exe
add_executable(resptest ../ut.cxx rres.cxx rfdr.cxx rmor.cxx
  rres_lda.cxx rres_gga.cxx rres_hgga.cxx rpp.cxx upp.cxx gpp.cxx
  misc.cxx)

target_compile_definitions(resptest PUBLIC BOOST_TEST_MODULE=RESPONSE)
target_include_directories(resptest PUBLIC ${RESP_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(resptest PUBLIC ${CQEX_LINK})

if(CQEX_DEP)
  add_dependencies(resptest ${CQEX_DEP})
endif()

# Generate directories
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/rresp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/rresp_ks)
#file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/uresp)
#file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/gresp)

file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/rresp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/rresp_ks)
#file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/uresp)
#file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/gresp)

file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/rmor)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/rmor_ks)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/rmor)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/rmor_ks)

file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/rpp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/upp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/gpp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/rpp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/upp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/gpp)

file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/misc)

# Add the Tests


add_test( RHF_RESIDUE   resptest --report_level=detailed --run_test=RHF_RESP/RHF_RESIDUE   )
#add_test( UHF_RESIDUE   resptest --report_level=detailed --run_test=UHF_RESP/UHF_RESIDUE   )
#add_test( GHF_RESIDUE   resptest --report_level=detailed --run_test=GHF_RESP/GHF_RESIDUE   )
#add_test( X2CHF_RESIDUE   resptest --report_level=detailed --run_test=X2CHF_RESP/X2CHF_RESIDUE   )
add_test( RKS_RESIDUE   resptest --report_level=detailed --run_test=RKS_RESP/RKS_RESIDUE   )
#add_test( UKS_RESIDUE   resptest --report_level=detailed --run_test=UKS_RESP/UKS_RESIDUE   )
#add_test( GKS_RESIDUE   resptest --report_level=detailed --run_test=GKS_RESP/GKS_RESIDUE   )
#add_test( X2CKS_RESIDUE   resptest --report_level=detailed --run_test=X2CKS_RESP/X2CKS_RESIDUE   )


add_test( RHF_PP_RESIDUE   resptest --report_level=detailed --run_test=RHF_RESP/RHF_PP_RESIDUE   )
add_test( UHF_PP_RESIDUE   resptest --report_level=detailed --run_test=UHF_RESP/UHF_PP_RESIDUE   )
#add_test( GHF_PP_RESIDUE   resptest --report_level=detailed --run_test=GHF_RESP/GHF_PP_RESIDUE   )
#add_test( X2CHF_PP_RESIDUE   resptest --report_level=detailed --run_test=X2CHF_RESP/X2CHF_PP_RESIDUE   )
#add_test( RKS_PP_RESIDUE   resptest --report_level=detailed --run_test=RKS_RESP/RKS_PP_RESIDUE   )
#add_test( UKS_PP_RESIDUE   resptest --report_level=detailed --run_test=UKS_RESP/UKS_PP_RESIDUE   )
#add_test( GKS_PP_RESIDUE   resptest --report_level=detailed --run_test=GKS_RESP/GKS_PP_RESIDUE   )
#add_test( X2CKS_PP_RESIDUE   resptest --report_level=detailed --run_test=X2CKS_RESP/X2CKS_PP_RESIDUE   )


add_test( RHF_FDR   resptest --report_level=detailed --run_test=RHF_RESP/RHF_FDR   )
#add_test( UHF_FDR   resptest --report_level=detailed --run_test=UHF_RESP/UHF_FDR   )
#add_test( GHF_FDR   resptest --report_level=detailed --run_test=GHF_RESP/GHF_FDR   )
#add_test( X2CHF_FDR   resptest --report_level=detailed --run_test=X2CHF_RESP/X2CHF_FDR   )
#add_test( RKS_FDR   resptest --report_level=detailed --run_test=RKS_RESP/RKS_FDR   )
#add_test( UKS_FDR   resptest --report_level=detailed --run_test=UKS_RESP/UKS_FDR   )
#add_test( GKS_FDR   resptest --report_level=detailed --run_test=GKS_RESP/GKS_FDR   )
#add_test( X2CKS_FDR   resptest --report_level=detailed --run_test=X2CKS_RESP/X2CKS_FDR   )


add_test( RHF_DFDR   resptest --report_level=detailed --run_test=RHF_RESP/RHF_DFDR   )
#add_test( UHF_DFDR   resptest --report_level=detailed --run_test=UHF_RESP/UHF_DFDR   )
#add_test( GHF_DFDR   resptest --report_level=detailed --run_test=GHF_RESP/GHF_DFDR   )
#add_test( X2CHF_DFDR   resptest --report_level=detailed --run_test=X2CHF_RESP/X2CHF_DFDR   )
#add_test( RKS_DFDR   resptest --report_level=detailed --run_test=RKS_RESP/RKS_DFDR   )
#add_test( UKS_DFDR   resptest --report_level=detailed --run_test=UKS_RESP/UKS_DFDR   )
#add_test( GKS_DFDR   resptest --report_level=detailed --run_test=GKS_RESP/GKS_DFDR   )
#add_test( X2CKS_DFDR   resptest --report_level=detailed --run_test=X2CKS_RESP/X2CKS_DFDR   )

add_test( RHF_MOR   resptest --report_level=detailed --run_test=RHF_MOR   )
#add_test( UHF_MOR   resptest --report_level=detailed --run_test=UHF_MOR   )
#add_test( GHF_MOR   resptest --report_level=detailed --run_test=GHF_MOR   )
#add_test( X2CHF_MOR   resptest --report_level=detailed --run_test=X2CHF_MOR   )
#add_test( RKS_MOR   resptest --report_level=detailed --run_test=RKS_MOR   )
#add_test( UKS_MOR   resptest --report_level=detailed --run_test=UKS_MOR   )
#add_test( GKS_MOR   resptest --report_level=detailed --run_test=GKS_MOR   )
#add_test( X2CKS_MOR   resptest --report_level=detailed --run_test=X2CKS_MOR   )

add_test ( MISC_RESP resptest --report_level=detailed --run_test=MISC_RESP )


if( CQ_ENABLE_MPI )

  add_test( NAME RHF_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RHF_RESP/RHF_RESIDUE" )
 #add_test( NAME UHF_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UHF_RESP/UHF_RESIDUE" )
 #add_test( NAME GHF_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GHF_RESP/GHF_RESIDUE" )
 #add_test( NAME X2CHF_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CHF_RESP/X2CHF_RESIDUE" )
  add_test( NAME RKS_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RKS_RESP/RKS_RESIDUE" )
 #add_test( NAME UKS_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UKS_RESP/UKS_RESIDUE" )
 #add_test( NAME GKS_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GKS_RESP/GKS_RESIDUE" )
 #add_test( NAME X2CKS_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CKS_RESP/X2CKS_RESIDUE" )


  add_test( NAME RHF_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RHF_RESP/RHF_FDR" )
 #add_test( NAME UHF_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UHF_RESP/UHF_FDR" )
 #add_test( NAME GHF_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GHF_RESP/GHF_FDR" )
 #add_test( NAME X2CHF_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CHF_RESP/X2CHF_FDR" )
 #add_test( NAME RKS_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RKS_RESP/RKS_FDR" )
 #add_test( NAME UKS_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UKS_RESP/UKS_FDR" )
 #add_test( NAME GKS_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GKS_RESP/GKS_FDR" )
 #add_test( NAME X2CKS_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CKS_RESP/X2CKS_FDR" )

  add_test( NAME RHF_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RHF_RESP/RHF_DFDR" )
 #add_test( NAME UHF_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UHF_RESP/UHF_DFDR" )
 #add_test( NAME GHF_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GHF_RESP/GHF_DFDR" )
 #add_test( NAME X2CHF_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CHF_RESP/X2CHF_DFDR" )
 #add_test( NAME RKS_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RKS_RESP/RKS_DFDR" )
 #add_test( NAME UKS_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UKS_RESP/UKS_DFDR" )
 #add_test( NAME GKS_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GKS_RESP/GKS_DFDR" )
 #add_test( NAME X2CKS_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CKS_RESP/X2CKS_DFDR" )



  add_test( NAME RHF_MOR_MPI      COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RHF_MOR"  )
 #add_test( NAME UHF_MOR_MPI      COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UHF_MOR"  )
 #add_test( NAME GHF_MOR_MPI      COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GHF_MOR"  )
 #add_test( NAME X2CHF_MOR_MPI      COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CHF_MOR"  )

endif()
