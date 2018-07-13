/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2018 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */

#include "resp.hpp"

#define CQRESTEST_IMPL(TNAME, JTYPE, IN, REF) \
BOOST_FIXTURE_TEST_CASE( TNAME, JTYPE ) { CQRESTEST( true, IN, REF ); }


BOOST_AUTO_TEST_SUITE( RHF_RESP )


BOOST_AUTO_TEST_SUITE( RHF_RESIDUE )



// FULL DIMENSIONAL TESTS

// Water 6-31G(d) TDHF (RESIDUE)
CQRESTEST_IMPL( Water_631Gd_RESIDUE, SerialJob, 
    "resp/serial/rresp/water_6-31Gd_rhf_residue",
    "water_6-31Gd_rhf_residue.bin.ref" )

#ifndef _CQ_GENERATE_TESTS

  // Water 6-31G(d) TDHF (RESIDUE, GPLHR)
  CQRESTEST_IMPL( Water_631Gd_RESIDUE_GPLHR, SerialJob, 
      "resp/serial/rresp/water_6-31Gd_rhf_residue_gplhr",
      "water_6-31Gd_rhf_residue.bin.ref" )
  
  // Water 6-31G(d) TDHF (RESIDUE, GPLHR + DIRECT)
  CQRESTEST_IMPL( Water_631Gd_RESIDUE_GPLHR_DIRECT, SerialJob, 
      "resp/serial/rresp/water_6-31Gd_rhf_residue_gplhr_direct",
      "water_6-31Gd_rhf_residue.bin.ref" )

#endif


#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) TDHF (RESIDUE)
  CQRESTEST_IMPL( PAR_Water_631Gd_RESIDUE, ParallelJob, 
      "resp/parallel/rresp/water_6-31Gd_rhf_residue",
      "water_6-31Gd_rhf_residue.bin.ref") 

  // SMP Water 6-31G(d) TDHF (RESIDUE, GPLHR)
  CQRESTEST_IMPL( PAR_Water_631Gd_RESIDUE_GPLHR, ParallelJob, 
      "resp/parallel/rresp/water_6-31Gd_rhf_residue_gplhr",
      "water_6-31Gd_rhf_residue.bin.ref" )
  
  // SMP Water 6-31G(d) TDHF (RESIDUE, GPLHR + DIRECT)
  CQRESTEST_IMPL( PAR_Water_631Gd_RESIDUE_GPLHR_DIRECT, ParallelJob, 
      "resp/parallel/rresp/water_6-31Gd_rhf_residue_gplhr_direct",
      "water_6-31Gd_rhf_residue.bin.ref" )

#endif


#ifndef _CQ_GENERATE_TESTS

  // A+B / A-B TESTS
    
  // Water 6-31G(d) TDHF (RESIDUE) A+B / A-B
  CQRESTEST_IMPL( Water_631Gd_RESIDUE_APB_AMB, SerialJob, 
      "resp/serial/rresp/water_6-31Gd_rhf_residue_apb_amb",
      "water_6-31Gd_rhf_residue.bin.ref") 
  

  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) TDHF (RESIDUE) A+B / A-B
    CQRESTEST_IMPL( PAR_Water_631Gd_RESIDUE_APB_AMB, ParallelJob, 
        "resp/parallel/rresp/water_6-31Gd_rhf_residue_apb_amb",
        "water_6-31Gd_rhf_residue.bin.ref") 
  

  #endif


  // REDUCED TESTS
    
  // Water 6-31G(d) TDHF (RESIDUE) REDUCED
  CQRESTEST_IMPL( Water_631Gd_RESIDUE_RED, SerialJob, 
      "resp/serial/rresp/water_6-31Gd_rhf_residue_reduced",
      "water_6-31Gd_rhf_residue.bin.ref") 
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) TDHF (RESIDUE) REDUCED
    CQRESTEST_IMPL( PAR_Water_631Gd_RESIDUE_RED, ParallelJob, 
        "resp/parallel/rresp/water_6-31Gd_rhf_residue_reduced",
        "water_6-31Gd_rhf_residue.bin.ref") 
  
  #endif


#endif


#if defined(CQ_ENABLE_MPI) && !defined(_CQ_GENERATE_TESTS)

  // Distributed-From-Root Full Dimensional Tests

  // Water 6-31G(d) TDHF (RESIDUE, GPLHR)
  CQRESTEST_IMPL( Water_631Gd_RESIDUE_GPLHR_DISTMATFROMROOT, SerialJob, 
      "resp/serial/rresp/water_6-31Gd_rhf_residue_gplhr_distmatfromroot",
      "water_6-31Gd_rhf_residue.bin.ref" )
  
  // SMP Water 6-31G(d) TDHF (RESIDUE, GPLHR)
  CQRESTEST_IMPL( PAR_Water_631Gd_RESIDUE_GPLHR_DISTMATFROMROOT, ParallelJob, 
      "resp/parallel/rresp/water_6-31Gd_rhf_residue_gplhr_distmatfromroot",
      "water_6-31Gd_rhf_residue.bin.ref" )





  // Build-Distributed Full Dimensional Tests

  // Water 6-31G(d) TDHF (RESIDUE, GPLHR)
  CQRESTEST_IMPL( Water_631Gd_RESIDUE_GPLHR_BUILDDIST, SerialJob, 
      "resp/serial/rresp/water_6-31Gd_rhf_residue_gplhr_builddist",
      "water_6-31Gd_rhf_residue.bin.ref" )
  
  // SMP Water 6-31G(d) TDHF (RESIDUE, GPLHR)
  CQRESTEST_IMPL( PAR_Water_631Gd_RESIDUE_GPLHR_BUILDDIST, ParallelJob, 
      "resp/parallel/rresp/water_6-31Gd_rhf_residue_gplhr_builddist",
      "water_6-31Gd_rhf_residue.bin.ref" )

#endif

// End RHF_RESIDUE suite
BOOST_AUTO_TEST_SUITE_END()

// End RHF_RESP suite
BOOST_AUTO_TEST_SUITE_END()


