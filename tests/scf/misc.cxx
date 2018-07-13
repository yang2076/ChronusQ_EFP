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

BOOST_AUTO_TEST_SUITE( MISC_SCF )

// Water 6-31G(d) {0., 0.01, 0.} Electric Field test
BOOST_FIXTURE_TEST_CASE( Water_631Gd_ed_0_0pt01_0, SerialJob ) {

  CQSCFTEST( "scf/serial/rhf/water_6-31Gd_ed_0_0.01_0", 
    "water_6-31Gd_ed_0_0.01_0.bin.ref" );
 
};

// Water 6-31G(d) {0.01, 0., 0.} Electric Field test
BOOST_FIXTURE_TEST_CASE( Water_631Gd_ed_0pt01_0_0, SerialJob ) {

  CQSCFTEST( "scf/serial/rhf/water_6-31Gd_ed_0.01_0_0", 
    "water_6-31Gd_ed_0.01_0_0.bin.ref" );
 
};

// Water 6-31G(d) {0., 0., 0.01} Electric Field test
BOOST_FIXTURE_TEST_CASE( Water_631Gd_ed_0_0_0pt01, SerialJob ) {

  CQSCFTEST( "scf/serial/rhf/water_6-31Gd_ed_0_0_0.01", 
    "water_6-31Gd_ed_0_0_0.01.bin.ref" );
 
};

// O2 Minimal basis
BOOST_FIXTURE_TEST_CASE( O2_STO3G, SerialJob ) {

  CQSCFTEST( "scf/serial/uhf/oxygen_sto-3g", "oxygen_sto-3g.bin.ref" );

};

BOOST_AUTO_TEST_SUITE_END()
