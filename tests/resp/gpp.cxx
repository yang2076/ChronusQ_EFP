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
 *    E-Mail: xsli@@uw.edu
 *  
 */

#include "resp.hpp"

#define CQRESTEST_IMPL(TNAME, JTYPE, IN, REF) \
BOOST_FIXTURE_TEST_CASE( TNAME, JTYPE ) { CQRESTEST( false, IN, REF ); }


BOOST_AUTO_TEST_SUITE( GHF_RESP )


BOOST_AUTO_TEST_SUITE( GHF_PP_RESIDUE )




// Li cc-pVDZ PP-RPA-HF (GHF) (RESIDUE)
CQRESTEST_IMPL( Li_ccpVDZ_ppRPA_GHF_RESIDUE, SerialJob, 
    "resp/serial/gpp/Li_ccpVDZ_GHF_ppRPA",
    "Li_ccpVDZ_GHF_ppRPA.bin.ref" )


#ifdef _CQ_DO_PARTESTS

  // SMP Li cc-pVDZ PP-RPA-HF (GHF) (RESIDUE)
  CQRESTEST_IMPL( PAR_Li_ccpVDZ_ppRPA_GHF_RESIDUE, ParallelJob, 
      "resp/parallel/gpp/Li_ccpVDZ_GHF_ppRPA",
      "Li_ccpVDZ_GHF_ppRPA.bin.ref" )

#endif




// End GHF_PP_RESIDUE suite
BOOST_AUTO_TEST_SUITE_END()

// End GHF_RESP suite
BOOST_AUTO_TEST_SUITE_END()
