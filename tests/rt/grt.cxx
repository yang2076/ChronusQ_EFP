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

#include "rt.hpp"


BOOST_AUTO_TEST_SUITE( X2CHF_RT )


// Water 6-311+G(d,p) X2CDelta Spike (along Y)
BOOST_FIXTURE_TEST_CASE( Water_6311pGdp_X2C_Delta_Y, SerialJob ) {

  CQRTTEST( rt/serial/grt/water_6-311pGdp_x2c_delta_y,
    water_6-311pGdp_x2c_delta_y.bin.ref );

}

#ifdef _CQ_DO_PARTESTS

// SMP Water 6-311+G(d,p) X2CDelta Spike (along Y)
BOOST_FIXTURE_TEST_CASE( PAR_Water_6311pGdp_X2C_Delta_Y, ParallelJob ) {

  CQRTTEST( rt/parallel/grt/water_6-311pGdp_x2c_delta_y,
    water_6-311pGdp_x2c_delta_y.bin.ref );

}

#endif


// End X2C_RT suite
BOOST_AUTO_TEST_SUITE_END()

