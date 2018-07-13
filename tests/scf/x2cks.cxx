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

BOOST_AUTO_TEST_SUITE( X2CKS )

// X2C Water 6-311+G(d,p) B3LYP (Spherical) test
BOOST_FIXTURE_TEST_CASE( Water_6311pGdp_x2c_b3lyp_sph, SerialJob ) {

  CQSCFTEST( "scf/serial/x2c/water_6-311+Gdp_b3lyp_sph", 
    "water_6-311+Gdp_sph_x2c_b3lyp.bin.ref",1e-6 );
 
};

// Water 6-311+G(d,p) B3LYP (Cartesian) test
BOOST_FIXTURE_TEST_CASE( Water_6311pGdp_x2c_b3lyp_cart, SerialJob ) {

  CQSCFTEST( "scf/serial/x2c/water_6-311+Gdp_b3lyp_cart", 
    "water_6-311+Gdp_cart_x2c_b3lyp.bin.ref",1e-6 );
 
};



// Hg SAPPORO DZP DKH_2012 SP SLATER
BOOST_FIXTURE_TEST_CASE( Hg_SAP_DZP_DKH3_2012_SP_SLATER , SerialJob ) {

  CQSCFTEST( "scf/serial/x2c/hg_sap_dz_dkh3_2012_sp_slater", 
    "hg_sap_dz_dkh3_2012_sp_slater.bin.ref",1e-6 );
 
};

// Zn SAPPORO DZP DKH_2012 SP SLATER
BOOST_FIXTURE_TEST_CASE( Zn_SAP_DZP_DKH3_2012_SP_SLATER , SerialJob ) {

  CQSCFTEST( "scf/serial/x2c/zn_sap_dz_dkh3_2012_sp_slater", 
    "zn_sap_dz_dkh3_2012_sp_slater.bin.ref",1e-6 );
 
};

// Cd SAPPORO DZP DKH_2012 SP SLATER
BOOST_FIXTURE_TEST_CASE( Cd_SAP_DZP_DKH3_2012_SP_SLATER , SerialJob ) {

  CQSCFTEST( "scf/serial/x2c/cd_sap_dz_dkh3_2012_sp_slater", 
    "cd_sap_dz_dkh3_2012_sp_slater.bin.ref",1e-6 );
 
};



// Hg SAPPORO DZP DKH_2012 SP B3LYP
BOOST_FIXTURE_TEST_CASE( Hg_SAP_DZP_DKH3_2012_SP_B3LYP , SerialJob ) {

  CQSCFTEST( "scf/serial/x2c/hg_sap_dz_dkh3_2012_sp_b3lyp", 
    "hg_sap_dz_dkh3_2012_sp_b3lyp.bin.ref",1e-6 );
 
};

// Zn SAPPORO DZP DKH_2012 SP B3LYP
BOOST_FIXTURE_TEST_CASE( Zn_SAP_DZP_DKH3_2012_SP_B3LYP , SerialJob ) {

  CQSCFTEST( "scf/serial/x2c/zn_sap_dz_dkh3_2012_sp_b3lyp", 
    "zn_sap_dz_dkh3_2012_sp_b3lyp.bin.ref",1e-6 );
 
};

// Cd SAPPORO DZP DKH_2012 SP B3LYP
BOOST_FIXTURE_TEST_CASE( Cd_SAP_DZP_DKH3_2012_SP_B3LYP , SerialJob ) {

  CQSCFTEST( "scf/serial/x2c/cd_sap_dz_dkh3_2012_sp_b3lyp", 
    "cd_sap_dz_dkh3_2012_sp_b3lyp.bin.ref",1e-6 );
 
};

#ifdef _CQ_DO_PARTESTS

// SMP X2C Water 6-311+G(d,p) B3LYP (Spherical) test
BOOST_FIXTURE_TEST_CASE( PAR_Water_6311pGdp_x2c_b3lyp_sph, ParallelJob ) {

  CQSCFTEST( "scf/parallel/x2c/water_6-311+Gdp_b3lyp_sph", 
    "water_6-311+Gdp_sph_x2c_b3lyp.bin.ref",1e-6 );
 
};


/*
// SMP Hg SAPPORO DZP DKH_2012 SP SLATER
BOOST_FIXTURE_TEST_CASE( PAR_Hg_SAP_DZP_DKH3_2012_SP_SLATER , ParallelJob ) {

  CQSCFTEST( scf/parallel/x2c/hg_sap_dz_dkh3_2012_sp_slater, 
    hg_sap_dz_dkh3_2012_sp_slater.bin.ref );
 
};

// SMP Hg SAPPORO DZP DKH_2012 SP B3LYP
BOOST_FIXTURE_TEST_CASE( PAR_Hg_SAP_DZP_DKH3_2012_SP_B3LYP , ParallelJob ) {

  CQSCFTEST( scf/parallel/x2c/hg_sap_dz_dkh3_2012_sp_b3lyp, 
    hg_sap_dz_dkh3_2012_sp_b3lyp.bin.ref );
 
};
*/

#endif

BOOST_AUTO_TEST_SUITE_END()



