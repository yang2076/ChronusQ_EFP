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

BOOST_AUTO_TEST_SUITE( X2CHF )

// Water 6-311+G(d,p) (Spherical) test
BOOST_FIXTURE_TEST_CASE( Water_6311pGdp_sph, SerialJob ) {

  CQSCFTEST( scf/serial/x2c/water_6-311+Gdp_sph, 
    water_6-311+Gdp_sph_x2c.bin.ref );
 
};

// Water 6-311+G(d,p) (Cartesian) test
BOOST_FIXTURE_TEST_CASE( Water_6311pGdp_cart, SerialJob ) {

  CQSCFTEST( scf/serial/x2c/water_6-311+Gdp_cart, 
    water_6-311+Gdp_cart_x2c.bin.ref );
 
};


#ifdef _CQ_DO_PARTESTS

// SMP Water 6-311+G(d,p) (Spherical) test
BOOST_FIXTURE_TEST_CASE( PAR_Water_6311pGdp_sph, ParallelJob ) {

  CQSCFTEST( scf/parallel/x2c/water_6-311+Gdp_sph, 
    water_6-311+Gdp_sph_x2c.bin.ref );
 
};


#endif

BOOST_AUTO_TEST_SUITE_END()



