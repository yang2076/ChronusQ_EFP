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


BOOST_AUTO_TEST_SUITE( UHF )

// Li 6-31G(d) test
BOOST_FIXTURE_TEST_CASE( Li_631Gd, SerialJob ) {

  CQSCFTEST( scf/serial/uhf/li_6-31Gd, li_6-31Gd.bin.ref );

};

#ifdef _CQ_DO_PARTESTS

// SMP Li 6-31G(d) test
BOOST_FIXTURE_TEST_CASE( PAR_Li_631Gd, ParallelJob ) {

  CQSCFTEST( scf/parallel/uhf/li_6-31Gd, li_6-31Gd.bin.ref );

};

#endif


BOOST_AUTO_TEST_SUITE_END()



