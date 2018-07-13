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

#ifndef _CQ_GENERATE_TESTS

#include "resp.hpp"

#define CQMORTEST_IMPL(TNAME, JTYPE, HER, AHER, IN, REF) \
BOOST_FIXTURE_TEST_CASE( TNAME, JTYPE ){ CQMORTEST(HER,AHER,IN,REF); }

BOOST_AUTO_TEST_SUITE( RHF_MOR )


// Full Dimensional Tests

// Water 6-31G(d) MOR
CQMORTEST_IMPL( Water_631Gd_MOR, SerialJob, true, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor",
    "water_6-31Gd_rhf_interp.bin.ref" )


// Water 6-31G(d) MOR (GMRES)
CQMORTEST_IMPL( Water_631Gd_MOR_GMRES, SerialJob, true, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_gmres",
    "water_6-31Gd_rhf_interp.bin.ref" )

// Water 6-31G(d) MOR (GMRES + DIRECT)
CQMORTEST_IMPL( Water_631Gd_MOR_GMRES_DIRECT, SerialJob, true, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_gmres_direct",
    "water_6-31Gd_rhf_interp.bin.ref" )

#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) MOR
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR, ParallelJob, true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_GMRES, ParallelJob, true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_gmres",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR (GMRES + DIRECT)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_GMRES_DIRECT, ParallelJob, true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_gmres_direct",
      "water_6-31Gd_rhf_interp.bin.ref" )

#endif





// A+B / A-B TESTS
   
// Water 6-31G(d) MOR A+B / A-B
CQMORTEST_IMPL( Water_631Gd_MOR_APB_AMB, SerialJob, true, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_apb_amb",
    "water_6-31Gd_rhf_interp.bin.ref" )


// Water 6-31G(d) MOR A+B / A-B (GMRES)
CQMORTEST_IMPL( Water_631Gd_MOR_APB_AMB_GMRES, SerialJob, true, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_apb_amb_gmres",
    "water_6-31Gd_rhf_interp.bin.ref" )

#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) MOR A+B /A-B
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_APB_AMB, ParallelJob, true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_apb_amb",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR A+B /A-B (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_APB_AMB_GMRES, ParallelJob, 
      true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_apb_amb_gmres",
      "water_6-31Gd_rhf_interp.bin.ref" )

#endif



// REDUCED HERMETIAN TESTS 
   
// Water 6-31G(d) MOR REDUCED HERMETIAN 
CQMORTEST_IMPL( Water_631Gd_MOR_RED_HER, SerialJob, true, false, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_her",
    "water_6-31Gd_rhf_interp.bin.ref" )

// Water 6-31G(d) MOR REDUCED HERMETIAN (GMRES)
CQMORTEST_IMPL( Water_631Gd_MOR_RED_HER_GMRES, SerialJob, true, false, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_her_gmres",
    "water_6-31Gd_rhf_interp.bin.ref" )

#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) MOR REDUCED HERMETIAN
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_HER, ParallelJob, true, false, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_her",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR REDUCED HERMETIAN (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_HER_GMRES, ParallelJob, 
      true, false, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_her_gmres",
      "water_6-31Gd_rhf_interp.bin.ref" )

#endif


#if 0

// REDUCED ANTI-HERMETIAN TESTS 
   
// Water 6-31G(d) MOR REDUCED ANTI-HERMETIAN 
CQMORTEST_IMPL( Water_631Gd_MOR_RED_ANTIHER, SerialJob, false, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_antiher",
    "water_6-31Gd_rhf_interp.bin.ref" )

// Water 6-31G(d) MOR REDUCED ANTI-HERMETIAN (GMRES)
CQMORTEST_IMPL( Water_631Gd_MOR_RED_ANTIHER_GMRES, SerialJob, false, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_antiher_gmres",
    "water_6-31Gd_rhf_interp.bin.ref" )

#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) MOR REDUCED ANTI-HERMETIAN
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_ANTIHER, ParallelJob, 
      false, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_antiher",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR REDUCED ANTI-HERMETIAN (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_ANTIHER_GMRES, ParallelJob, 
      false, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_antiher_gmres",
      "water_6-31Gd_rhf_interp.bin.ref" )

#endif

#endif





#if defined(CQ_ENABLE_MPI)

  // Distributed Full Dimensional Tests

  // Water 6-31G(d) MOR
  CQMORTEST_IMPL( Water_631Gd_MOR_DISTMATFROMROOT, SerialJob, true, true, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_DISTMATFROMROOT, ParallelJob, 
      true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // Water 6-31G(d) MOR (GMRES)
  CQMORTEST_IMPL( Water_631Gd_MOR_GMRES_DISTMATFROMROOT, SerialJob, 
      true, true, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_GMRES_DISTMATFROMROOT, ParallelJob, 
      true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )




  // Distributed A+B / A-B Dimensional Tests
    
  // Water 6-31G(d) MOR A+B / A-B
  CQMORTEST_IMPL( Water_631Gd_MOR_APB_AMB_DISTMATFROMROOT, SerialJob, 
      true, true, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_apb_amb_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR A+B / A-B
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_APB_AMB_DISTMATFROMROOT, ParallelJob, 
      true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_apb_amb_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // Water 6-31G(d) MOR A+B / A-B (GMRES)
  CQMORTEST_IMPL( Water_631Gd_MOR_APB_AMB_GMRES_DISTMATFROMROOT, SerialJob, 
      true, true, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_apb_amb_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR A+B / A-B (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_APB_AMB_GMRES_DISTMATFROMROOT, 
      ParallelJob, 
      true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_apb_amb_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )



  // Distributed Reduced Hermetian Tests

  // Water 6-31G(d) MOR REDUCED HERMETIAN
  CQMORTEST_IMPL( Water_631Gd_MOR_RED_HER_DISTMATFROMROOT, SerialJob, 
      true, false, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_her_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR REDUCED HERMETIAN
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_HER_DISTMATFROMROOT, ParallelJob, 
      true, false, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_her_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // Water 6-31G(d) MOR REDUCED HERMETIAN (GMRES)
  CQMORTEST_IMPL( Water_631Gd_MOR_RED_HER_GMRES_DISTMATFROMROOT, SerialJob, 
      true, false, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_her_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR REDUCED HERMETIAN (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_HER_GMRES_DISTMATFROMROOT, 
      ParallelJob, 
      true, false, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_her_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

#if 0


  // Distributed Reduced Anti-Hermetian Tests

  // Water 6-31G(d) MOR REDUCED ANTIHERMETIAN
  CQMORTEST_IMPL( Water_631Gd_MOR_RED_ANTIHER_DISTMATFROMROOT, SerialJob, 
      false, true, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_antiher_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR REDUCED ANTIHERMETIAN
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_ANTIHER_DISTMATFROMROOT, 
      ParallelJob, false, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_antiher_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // Water 6-31G(d) MOR REDUCED ANTIHERMETIAN (GMRES)
  CQMORTEST_IMPL( Water_631Gd_MOR_RED_ANTIHER_GMRES_DISTMATFROMROOT, 
      SerialJob, 
      false, true, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_antiher_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR REDUCED ANTIHERMETIAN (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_ANTIHER_GMRES_DISTMATFROMROOT, 
      ParallelJob, false, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_antiher_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )


#endif

#endif


// End RHF_MOR suite
BOOST_AUTO_TEST_SUITE_END()

#endif
