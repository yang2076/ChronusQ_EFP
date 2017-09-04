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
  BOOST_CHECK_MESSAGE(std::abs(yDummy - xDummy) < 9e-10, "ENERGY TEST FAILED " << std::abs(yDummy - xDummy) );\
  \
  /* Check Multipoles */ \
  std::array<double,3> xDummy3, yDummy3; \
  refFile.readData("SCF/LEN_ELECTRIC_DIPOLE",&xDummy3[0]);\
  resFile.readData("SCF/LEN_ELECTRIC_DIPOLE",&yDummy3[0]);\
  for(auto i = 0; i < 3; i++)\
    BOOST_CHECK_MESSAGE(std::abs(yDummy3[i] - xDummy3[i]) < 9e-9, "DIPOLE TEST FAILED IXYZ = " << i << " " << std::abs(yDummy3[i] - xDummy3[i])  );\
  \
  std::array<std::array<double,3>,3> xDummy33, yDummy33;\
  refFile.readData("SCF/LEN_ELECTRIC_QUADRUPOLE",&xDummy33[0][0]);\
  resFile.readData("SCF/LEN_ELECTRIC_QUADRUPOLE",&yDummy33[0][0]);\
  for(auto i = 0; i < 3; i++)\
  for(auto j = 0; j < 3; j++)\
    BOOST_CHECK_MESSAGE(std::abs(yDummy33[i][j] - xDummy33[i][j]) < 9e-9, "QUADRUPOLE TEST FAILED IXYZ = " << i << " JXYZ = " << j << " " << std::abs(yDummy33[i][j] - xDummy33[i][j]) );\
  \
  std::array<std::array<std::array<double,3>,3>,3> xDummy333, yDummy333;\
  refFile.readData("SCF/LEN_ELECTRIC_OCTUPOLE",&xDummy333[0][0][0]);\
  resFile.readData("SCF/LEN_ELECTRIC_OCTUPOLE",&yDummy333[0][0][0]);\
  for(auto i = 0; i < 3; i++)\
  for(auto j = 0; j < 3; j++)\
  for(auto k = 0; k < 3; k++)\
    BOOST_CHECK_MESSAGE(std::abs(yDummy333[i][j][k] - xDummy333[i][j][k]) < 9e-9, "OCTUPOLE TEST FAILED IXYZ = " << i << " JXYZ = " << j << " KXYZ = " << k << " " << std::abs(yDummy333[i][j][k] - xDummy333[i][j][k]) );\
  \
  /* Check Spin */\
  refFile.readData("SCF/S_EXPECT",&xDummy3[0]);\
  resFile.readData("SCF/S_EXPECT",&yDummy3[0]);\
  for(auto i = 0; i < 3; i++)\
    BOOST_CHECK_MESSAGE(std::abs(yDummy3[i] - xDummy3[i]) < 9e-10, "<S> TEST FAILED IXYZ = " << i << " " << std::abs(yDummy3[i] - xDummy3[i]) );\
  \
  refFile.readData("SCF/S_SQUARED",&xDummy);\
  resFile.readData("SCF/S_SQUARED",&yDummy);\
  BOOST_CHECK_MESSAGE(std::abs(yDummy - xDummy) < 9e-10, "<S^2> TEST FAILED " << std::abs(yDummy - xDummy) );
  
  

#endif

#endif
