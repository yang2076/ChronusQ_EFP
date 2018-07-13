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

#include "scf.hpp"


BOOST_AUTO_TEST_SUITE( RKS )


// B3LYP / cc-pVTZ
BOOST_FIXTURE_TEST_CASE( Water_ccpVTZ_B3LYP, SerialJob ) {

  CQSCFTEST( "scf/serial/rks/water_cc-pVTZ_B3LYP", "water_cc-pVTZ_B3LYP.bin.ref" );

}

// BLYP / cc-pVTZ
BOOST_FIXTURE_TEST_CASE( Water_ccpVTZ_BLYP, SerialJob ) {

  CQSCFTEST( "scf/serial/rks/water_cc-pVTZ_BLYP", "water_cc-pVTZ_BLYP.bin.ref" );

}

// LSDA / cc-pVTZ
BOOST_FIXTURE_TEST_CASE( Water_ccpVTZ_LSDA, SerialJob ) {

  CQSCFTEST( "scf/serial/rks/water_cc-pVTZ_LSDA", "water_cc-pVTZ_LSDA.bin.ref" );

}

#ifdef _CQ_DO_PARTESTS

// SMP B3LYP / cc-pVTZ
BOOST_FIXTURE_TEST_CASE( PAR_Water_ccpVTZ_B3LYP, ParallelJob ) {

  CQSCFTEST( "scf/parallel/rks/water_cc-pVTZ_B3LYP", "water_cc-pVTZ_B3LYP.bin.ref" );

}

// SMP BLYP / cc-pVTZ
BOOST_FIXTURE_TEST_CASE( PAR_Water_ccpVTZ_BLYP, ParallelJob ) {

  CQSCFTEST( "scf/parallel/rks/water_cc-pVTZ_BLYP", "water_cc-pVTZ_BLYP.bin.ref" );

}

// SMP LSDA / cc-pVTZ
BOOST_FIXTURE_TEST_CASE( PAR_Water_ccpVTZ_LSDA, ParallelJob ) {

  CQSCFTEST( "scf/parallel/rks/water_cc-pVTZ_LSDA", "water_cc-pVTZ_LSDA.bin.ref" );

}

#endif


// End RKS suite
BOOST_AUTO_TEST_SUITE_END()


