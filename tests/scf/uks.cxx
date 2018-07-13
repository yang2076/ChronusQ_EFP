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


BOOST_AUTO_TEST_SUITE( UKS )


// B3LYP / 6-311pG**
BOOST_FIXTURE_TEST_CASE( Oxygen_6311pGss_B3LYP, SerialJob ) {

  CQSCFTEST( "scf/serial/uks/oxygen_6-311pG**_B3LYP", "oxygen_6-311pG**_B3LYP.bin.ref" );

}

// BLYP / 6-311pG**
BOOST_FIXTURE_TEST_CASE( Oxygen_6311pGss_BLYP, SerialJob ) {

  CQSCFTEST( "scf/serial/uks/oxygen_6-311pG**_BLYP", "oxygen_6-311pG**_BLYP.bin.ref" );

}

// LSDA / 6-311pG**
BOOST_FIXTURE_TEST_CASE( Oxygen_6311pGss_LSDA, SerialJob ) {

  CQSCFTEST( "scf/serial/uks/oxygen_6-311pG**_LSDA", "oxygen_6-311pG**_LSDA.bin.ref" );

}

#ifdef _CQ_DO_PARTESTS

// SMP B3LYP / 6-311pG**
BOOST_FIXTURE_TEST_CASE( PAR_Oxygen_6311pGss_B3LYP, ParallelJob ) {

  CQSCFTEST( "scf/parallel/uks/oxygen_6-311pG**_B3LYP", "oxygen_6-311pG**_B3LYP.bin.ref" );

}

// SMP BLYP / 6-311pG**
BOOST_FIXTURE_TEST_CASE( PAR_Oxygen_6311pGss_BLYP, ParallelJob ) {

  CQSCFTEST( "scf/parallel/uks/oxygen_6-311pG**_BLYP", "oxygen_6-311pG**_BLYP.bin.ref" );

}

// SMP LSDA / 6-311pG**
BOOST_FIXTURE_TEST_CASE( PAR_Oxygen_6311pGss_LSDA, ParallelJob ) {

  CQSCFTEST( "scf/parallel/uks/oxygen_6-311pG**_LSDA", "oxygen_6-311pG**_LSDA.bin.ref" );

}

#endif


// End UKS suite
BOOST_AUTO_TEST_SUITE_END()


