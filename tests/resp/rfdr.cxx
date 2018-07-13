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


#define CQFDRTEST_IMPL(TNAME, JTYPE, HER, AHER, IN, REF) \
BOOST_FIXTURE_TEST_CASE( TNAME, JTYPE ){ CQFDRTEST<double>(HER,AHER,IN,REF); }

#define CQDFDRTEST_IMPL(TNAME, JTYPE, HER, AHER, IN, REF) \
BOOST_FIXTURE_TEST_CASE( TNAME, JTYPE ){ CQFDRTEST<dcomplex>(HER,AHER,IN,REF); }

BOOST_AUTO_TEST_SUITE( RHF_RESP )


BOOST_AUTO_TEST_SUITE( RHF_FDR )

// FULL DIMENSIONAL TESTS
  
// Water 6-31G(d) FDR
CQFDRTEST_IMPL( Water_631Gd_FDR, SerialJob, true, true, 
    "resp/serial/rresp/water_6-31Gd_rhf_fdr",
    "water_6-31Gd_rhf_fdr.bin.ref" )

#ifndef _CQ_GENERATE_TESTS

  // Water 6-31G(d) FDR (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_GMRES, SerialJob, true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_gmres",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // Water 6-31G(d) FDR (GMRES + DIRECT)
  CQFDRTEST_IMPL( Water_631Gd_FDR_GMRES_DIRECT, SerialJob, true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_gmres_direct",
      "water_6-31Gd_rhf_fdr.bin.ref" )

#endif

#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) FDR
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR, ParallelJob, true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR (GMRES)
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_GMRES, ParallelJob, true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_gmres",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR (GMRES + DIRECT)
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_GMRES_DIRECT, ParallelJob, true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_gmres_direct",
      "water_6-31Gd_rhf_fdr.bin.ref" )

#endif


#ifndef _CQ_GENERATE_TESTS

  // A+B / A-B TESTS
  
  // Water 6-31G(d) FDR A+B / A-B
  CQFDRTEST_IMPL( Water_631Gd_FDR_APB_AMB, SerialJob, true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_apb_amb",
      "water_6-31Gd_rhf_fdr.bin.ref" )


  // Water 6-31G(d) FDR A+B / A-B (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_APB_AMB_GMRES, SerialJob, true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_apb_amb_gmres",
      "water_6-31Gd_rhf_fdr.bin.ref" )
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) FDR A+B / A-B
    CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_APB_AMB, ParallelJob, true, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_fdr_apb_amb",
        "water_6-31Gd_rhf_fdr.bin.ref" )


    // SMP Water 6-31G(d) FDR A+B / A-B (GMRES)
    CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_APB_AMB_GMRES, ParallelJob, true, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_fdr_apb_amb_gmres",
        "water_6-31Gd_rhf_fdr.bin.ref" )
  
  #endif



  // REDUCED HERMETIAN TESTS
  
  // Water 6-31G(d) FDR REDUCED HERMETIAN
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_HER, SerialJob, true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_her",
      "water_6-31Gd_rhf_fdr.bin.ref" )


  // Water 6-31G(d) FDR REDUCED HERMETIAN (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_HER_GMRES, SerialJob, true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_her_gmres",
      "water_6-31Gd_rhf_fdr.bin.ref" )
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) FDR REDUCED HERMETIAN
    CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_HER, ParallelJob, true, false, 
        "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_her",
        "water_6-31Gd_rhf_fdr.bin.ref" )

    // SMP Water 6-31G(d) FDR REDUCED HERMETIAN (GMRES)
    CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_HER_GMRES, ParallelJob, 
        true, false, 
        "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_her_gmres",
        "water_6-31Gd_rhf_fdr.bin.ref" )
  
  #endif




  // REDUCED ANTI-HERMETIAN TESTS
  
  // Water 6-31G(d) FDR REDUCED ANTI-HERMETIAN
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_ANTIHER, SerialJob, false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_antiher",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // Water 6-31G(d) FDR REDUCED ANTI-HERMETIAN (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_ANTIHER_GMRES, SerialJob, false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_antiher_gmres",
      "water_6-31Gd_rhf_fdr.bin.ref" )
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) FDR REDUCED ANTI-HERMETIAN
    CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_ANTIHER, ParallelJob, false, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_antiher",
        "water_6-31Gd_rhf_fdr.bin.ref" )

    // SMP Water 6-31G(d) FDR REDUCED ANTI-HERMETIAN (GMRES)
    CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_ANTIHER_GMRES, ParallelJob, 
        false, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_antiher_gmres",
        "water_6-31Gd_rhf_fdr.bin.ref" )
  
  #endif

#endif


#if defined(CQ_ENABLE_MPI) && !defined(_CQ_GENERATE_TESTS)

  // Distributed Full Dimensional Tests

  // Water 6-31G(d) FDR
  CQFDRTEST_IMPL( Water_631Gd_FDR_DISTMATFROMROOT, SerialJob, true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_DISTMATFROMROOT, ParallelJob, true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // Water 6-31G(d) FDR (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_GMRES_DISTMATFROMROOT, SerialJob, true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR (GMRES)
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_GMRES_DISTMATFROMROOT, ParallelJob, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )




  // Distributed A+B / A-B Dimensional Tests
    
  // Water 6-31G(d) FDR A+B / A-B
  CQFDRTEST_IMPL( Water_631Gd_FDR_APB_AMB_DISTMATFROMROOT, SerialJob, 
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_apb_amb_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR A+B / A-B
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_APB_AMB_DISTMATFROMROOT, ParallelJob, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_apb_amb_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // Water 6-31G(d) FDR A+B / A-B (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_APB_AMB_GMRES_DISTMATFROMROOT, SerialJob, 
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_apb_amb_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR A+B / A-B (GMRES)
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_APB_AMB_GMRES_DISTMATFROMROOT, 
      ParallelJob, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_apb_amb_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )



  // Distributed Reduced Hermetian Tests

  // Water 6-31G(d) FDR REDUCED HERMETIAN
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_HER_DISTMATFROMROOT, SerialJob, 
      true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_her_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR REDUCED HERMETIAN
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_HER_DISTMATFROMROOT, ParallelJob, 
      true, false, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_her_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // Water 6-31G(d) FDR REDUCED HERMETIAN (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_HER_GMRES_DISTMATFROMROOT, SerialJob, 
      true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_her_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR REDUCED HERMETIAN (GMRES)
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_HER_GMRES_DISTMATFROMROOT, 
      ParallelJob, 
      true, false, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_her_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )



  // Distributed Reduced Anti-Hermetian Tests

  // Water 6-31G(d) FDR REDUCED ANTIHERMETIAN
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_ANTIHER_DISTMATFROMROOT, SerialJob, 
      false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_antiher_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR REDUCED ANTIHERMETIAN
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_ANTIHER_DISTMATFROMROOT, 
      ParallelJob, false, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_antiher_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // Water 6-31G(d) FDR REDUCED ANTIHERMETIAN (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_ANTIHER_GMRES_DISTMATFROMROOT, SerialJob, 
      false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_antiher_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR REDUCED ANTIHERMETIAN (GMRES)
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_ANTIHER_GMRES_DISTMATFROMROOT, 
      ParallelJob, false, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_antiher_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

#endif



// End RHF_FDR suite
BOOST_AUTO_TEST_SUITE_END()

























BOOST_AUTO_TEST_SUITE( RHF_DFDR )


// Full Dimensional Tests

// Water 6-31G(d) DFDR
CQDFDRTEST_IMPL( Water_631Gd_DFDR, SerialJob, true, true, 
    "resp/serial/rresp/water_6-31Gd_rhf_dfdr",
    "water_6-31Gd_rhf_dfdr.bin.ref" )

#ifndef _CQ_GENERATE_TESTS

  // Water 6-31G(d) DFDR (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_GMRES, SerialJob, true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_gmres",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR (GMRES + DIRECT)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_GMRES_DIRECT, SerialJob, true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_gmres_direct",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

#endif

#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) DFDR
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR, ParallelJob, true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_GMRES, ParallelJob, true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_gmres",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR (GMRES + DIRECT)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_GMRES_DIRECT, ParallelJob, true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_gmres_direct",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

#endif




#ifndef _CQ_GENERATE_TESTS

  // A+B / A-B TESTS
     
  // Water 6-31G(d) DFDR A+B / A-B
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_APB_AMB, SerialJob, true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_apb_amb",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR A+B / A-B (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_APB_AMB_GMRES, SerialJob, true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_apb_amb_gmres",
      "water_6-31Gd_rhf_dfdr.bin.ref" )
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) DFDR A+B /A-B
    CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_APB_AMB, ParallelJob, true, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_apb_amb",
        "water_6-31Gd_rhf_dfdr.bin.ref" )

    // SMP Water 6-31G(d) DFDR A+B /A-B (GMRES)
    CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_APB_AMB_GMRES, ParallelJob, 
        true, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_apb_amb_gmres",
        "water_6-31Gd_rhf_dfdr.bin.ref" )
  
  #endif


  // REDUCED HERMETIAN TESTS 
     
  // Water 6-31G(d) DFDR REDUCED HERMETIAN 
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_HER, SerialJob, true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_her",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR REDUCED HERMETIAN (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_HER_GMRES, SerialJob, true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_her_gmres",
      "water_6-31Gd_rhf_dfdr.bin.ref" )
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) DFDR REDUCED HERMETIAN
    CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_HER, ParallelJob, true, false, 
        "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_her",
        "water_6-31Gd_rhf_dfdr.bin.ref" )

    // SMP Water 6-31G(d) DFDR REDUCED HERMETIAN (GMRES)
    CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_HER_GMRES, ParallelJob, 
        true, false, 
        "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_her_gmres",
        "water_6-31Gd_rhf_dfdr.bin.ref" )
  
  #endif



  // REDUCED ANTI-HERMETIAN TESTS 
     
  // Water 6-31G(d) DFDR REDUCED ANTI-HERMETIAN 
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_ANTIHER, SerialJob, false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR REDUCED ANTI-HERMETIAN (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_ANTIHER_GMRES, SerialJob, false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_gmres",
      "water_6-31Gd_rhf_dfdr.bin.ref" )
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) DFDR REDUCED ANTI-HERMETIAN
    CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_ANTIHER, ParallelJob, 
        false, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher",
        "water_6-31Gd_rhf_dfdr.bin.ref" )

    // SMP Water 6-31G(d) DFDR REDUCED ANTI-HERMETIAN (GMRES)
    CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_ANTIHER_GMRES, ParallelJob, 
        false, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_gmres",
        "water_6-31Gd_rhf_dfdr.bin.ref" )
  
  #endif


#endif



#if defined(CQ_ENABLE_MPI) && !defined(_CQ_GENERATE_TESTS)

  // Distributed-From-Root Full Dimensional Tests

  // Water 6-31G(d) DFDR
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_DISTMATFROMROOT, SerialJob, true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_DISTMATFROMROOT, ParallelJob, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_GMRES_DISTMATFROMROOT, SerialJob, 
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_GMRES_DISTMATFROMROOT, ParallelJob, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )




  // Distributed-From-Root A+B / A-B Dimensional Tests
    
  // Water 6-31G(d) DFDR A+B / A-B
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_APB_AMB_DISTMATFROMROOT, SerialJob, 
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_apb_amb_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR A+B / A-B
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_APB_AMB_DISTMATFROMROOT, ParallelJob, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_apb_amb_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR A+B / A-B (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_APB_AMB_GMRES_DISTMATFROMROOT, SerialJob, 
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_apb_amb_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR A+B / A-B (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_APB_AMB_GMRES_DISTMATFROMROOT, 
      ParallelJob, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_apb_amb_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )



  // Distributed-From-Root Reduced Hermetian Tests

  // Water 6-31G(d) DFDR REDUCED HERMETIAN
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_HER_DISTMATFROMROOT, SerialJob, 
      true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_her_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED HERMETIAN
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_HER_DISTMATFROMROOT, ParallelJob, 
      true, false, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_her_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR REDUCED HERMETIAN (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_HER_GMRES_DISTMATFROMROOT, SerialJob, 
      true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_her_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED HERMETIAN (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_HER_GMRES_DISTMATFROMROOT, 
      ParallelJob, 
      true, false, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_her_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )



  // Distributed-From-Root Reduced Anti-Hermetian Tests

  // Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_ANTIHER_DISTMATFROMROOT, SerialJob, 
      false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_ANTIHER_DISTMATFROMROOT, 
      ParallelJob, false, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_ANTIHER_GMRES_DISTMATFROMROOT, 
      SerialJob, 
      false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_ANTIHER_GMRES_DISTMATFROMROOT, 
      ParallelJob, false, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )












  // Build-Distributed Full Dimensional Tests

  // Water 6-31G(d) DFDR
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_BUILDDIST, SerialJob, true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_BUILDDIST, ParallelJob, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_GMRES_BUILDDIST, SerialJob, 
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_GMRES_BUILDDIST, ParallelJob, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )




  // Build-Distributed A+B / A-B Dimensional Tests
    
  // Water 6-31G(d) DFDR A+B / A-B
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_APB_AMB_BUILDDIST, SerialJob, 
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_apb_amb_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR A+B / A-B
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_APB_AMB_BUILDDIST, ParallelJob, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_apb_amb_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR A+B / A-B (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_APB_AMB_GMRES_BUILDDIST, SerialJob, 
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_apb_amb_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR A+B / A-B (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_APB_AMB_GMRES_BUILDDIST, 
      ParallelJob, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_apb_amb_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )



  // Build-Distributed Reduced Hermetian Tests

  // Water 6-31G(d) DFDR REDUCED HERMETIAN
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_HER_BUILDDIST, SerialJob, 
      true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_her_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED HERMETIAN
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_HER_BUILDDIST, ParallelJob, 
      true, false, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_her_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR REDUCED HERMETIAN (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_HER_GMRES_BUILDDIST, SerialJob, 
      true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_her_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED HERMETIAN (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_HER_GMRES_BUILDDIST, 
      ParallelJob, 
      true, false, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_her_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )



  // Build-Distributed Reduced Anti-Hermetian Tests

  // Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_ANTIHER_BUILDDIST, SerialJob, 
      false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_ANTIHER_BUILDDIST, 
      ParallelJob, false, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_ANTIHER_GMRES_BUILDDIST, 
      SerialJob, 
      false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_ANTIHER_GMRES_BUILDDIST, 
      ParallelJob, false, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

#endif

// End RHF_DFDR suite
BOOST_AUTO_TEST_SUITE_END()




// End RHF_RESP suite
BOOST_AUTO_TEST_SUITE_END()
