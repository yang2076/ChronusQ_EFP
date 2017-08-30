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

#ifndef __INCLUDED_TESTS_SCF_HPP__
#define __INCLUDED_TESTS_SCF_HPP__

#include <ut.hpp>

#include <cxxapi/procedural.hpp>
#include <util/files.hpp>

// Directory containing reference files
#define SCF_TEST_REF TEST_ROOT "/scf/reference/"

using namespace ChronusQ;


#ifdef _CQ_GENERATE_TESTS

// Run CQ job and write reference files
// *** WARNING: This will overwrite existing reference files ***
// ***                       USE AOR                         ***
#define CQSCFTEST( in, ref ) \
  RunChronusQ(TEST_ROOT #in ".inp","STDOUT", \
    SCF_TEST_REF #ref,TEST_OUT #in ".scr");


#else

// HTG SCF test
#define CQSCFTEST( in, ref ) \
  RunChronusQ(TEST_ROOT #in ".inp","STDOUT", \
    TEST_OUT #in ".bin",TEST_OUT #in ".scr");\
  \
  SafeFile refFile(SCF_TEST_REF #ref,true);\
  SafeFile resFile(TEST_OUT #in ".bin",true);\
  \
  double xDummy, yDummy;\
  \
  /* Check Energy */ \
  refFile.readData("SCF/TOTAL_ENERGY",&xDummy);\
  resFile.readData("SCF/TOTAL_ENERGY",&yDummy);\
  BOOST_CHECK((yDummy - xDummy) < 1e-10);

#endif

#endif
