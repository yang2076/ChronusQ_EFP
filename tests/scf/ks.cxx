/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2017 Li Research Group (University of Washington)
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


BOOST_AUTO_TEST_SUITE( KS_KEYWORD )


// B3LYP
BOOST_FIXTURE_TEST_CASE( KEYWORD_B3LYP, SerialJob ) {

  CQSCFTEST( scf/serial/rks/water_sto-3g_B3LYP, water_sto-3g_B3LYP.bin.ref );

}

// BHANDH
BOOST_FIXTURE_TEST_CASE( KEYWORD_BHANDH, SerialJob ) {

  CQSCFTEST( scf/serial/rks/water_sto-3g_BHANDH, water_sto-3g_BHANDH.bin.ref );

}

// BHANDHLYP
BOOST_FIXTURE_TEST_CASE( KEYWORD_BHANDHLYP, SerialJob ) {

  CQSCFTEST( scf/serial/rks/water_sto-3g_BHANDHLYP, water_sto-3g_BHANDHLYP.bin.ref );

}

// BLYP
BOOST_FIXTURE_TEST_CASE( KEYWORD_BLYP, SerialJob ) {

  CQSCFTEST( scf/serial/rks/water_sto-3g_BLYP, water_sto-3g_BLYP.bin.ref );

}

// LSDA
BOOST_FIXTURE_TEST_CASE( KEYWORD_LSDA, SerialJob ) {

  CQSCFTEST( scf/serial/rks/water_sto-3g_LSDA, water_sto-3g_LSDA.bin.ref );

}

// PBE0
BOOST_FIXTURE_TEST_CASE( KEYWORD_PBE0, SerialJob ) {

  CQSCFTEST( scf/serial/rks/water_sto-3g_PBE0, water_sto-3g_PBE0.bin.ref );

}

// PBEXPBEC
BOOST_FIXTURE_TEST_CASE( KEYWORD_PBEXPBEC, SerialJob ) {

  CQSCFTEST( scf/serial/rks/water_sto-3g_PBEXPBEC, water_sto-3g_PBEXPBEC.bin.ref );

}

// SLATER
BOOST_FIXTURE_TEST_CASE( KEYWORD_SLATER, SerialJob ) {

  CQSCFTEST( scf/serial/rks/water_sto-3g_SLATER, water_sto-3g_SLATER.bin.ref );

}

// End KS_KEYWORD suite
BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE( KS_FUNC )

// B3LYP
BOOST_FIXTURE_TEST_CASE( KS_CART_B3LYP, SerialJob ) {

  CQSCFTEST( scf/serial/rks/water_cc-pVTZ_cart_B3LYP, water_cc-pVTZ_cart_B3LYP.bin.ref );

}

// End KS_FUNC suite
BOOST_AUTO_TEST_SUITE_END()


