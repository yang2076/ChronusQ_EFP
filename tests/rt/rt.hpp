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

#ifndef __INCLUDED_TESTS_RT_HPP__
#define __INCLUDED_TESTS_RT_HPP__

#include <ut.hpp>

#include <cxxapi/procedural.hpp>
#include <util/files.hpp>

// Directory containing reference files
#define RT_TEST_REF TEST_ROOT "/rt/reference/"

using namespace ChronusQ;


#ifdef _CQ_GENERATE_TESTS

// Run CQ job and write reference files
// *** WARNING: This will overwrite existing reference files ***
// ***                       USE AOR                         ***
#define CQRTTEST( in, ref ) \
  RunChronusQ(TEST_ROOT #in ".inp","STDOUT", \
    RT_TEST_REF #ref,TEST_OUT #in ".scr");


#else

// HTG RT test
#define CQRTTEST( in, ref ) \
  RunChronusQ(TEST_ROOT #in ".inp","STDOUT", \
    TEST_OUT #in ".bin",TEST_OUT #in ".scr");\
  \
  SafeFile refFile(RT_TEST_REF #ref,true);\
  SafeFile resFile(TEST_OUT #in ".bin",true);\
  \
  std::vector<double> xDummy, yDummy;\
  std::vector<std::array<double,3>> xDummy3, yDummy3; \
  \
  auto energyDim1 = resFile.getDims("/RT/ENERGY");\
  auto energyDim2 = refFile.getDims("/RT/ENERGY");\
  std::cerr << energyDim1[0] << ", " << energyDim2[0] << std::endl;\
  if( energyDim1.size() != 1 or energyDim2.size() != 1 or \
      energyDim1[0] != energyDim2[0] ) \
    BOOST_FAIL("Something went wrong in the file generation for energies");\
  \
  auto dipoleDim1 = resFile.getDims("/RT/LEN_ELEC_DIPOLE");\
  auto dipoleDim2 = refFile.getDims("/RT/LEN_ELEC_DIPOLE");\
  if( dipoleDim1.size() != 2 or dipoleDim2.size() != 2 or \
      dipoleDim1[0] != dipoleDim2[0] or dipoleDim1[1] != 3 or \
      dipoleDim2[1] != 3) \
    BOOST_FAIL("Something went wrong in the file generation for dipoles");\
  \
  xDummy.resize(energyDim1[0]); yDummy.resize(energyDim1[0]);\
  resFile.readData("/RT/ENERGY",&xDummy[0]);\
  refFile.readData("/RT/ENERGY",&yDummy[0]);\
  \
  for(auto i = 0; i < energyDim1[0]; i++) \
    BOOST_CHECK(std::abs(xDummy[i] - yDummy[i]) < 1e-9);\
  \
  xDummy3.resize(dipoleDim1[0]); yDummy3.resize(dipoleDim1[0]);\
  resFile.readData("/RT/LEN_ELEC_DIPOLE",&xDummy3[0][0]);\
  refFile.readData("/RT/LEN_ELEC_DIPOLE",&yDummy3[0][0]);\
  \
  for(auto i = 0; i < energyDim1[0]; i++) {\
    BOOST_CHECK(std::abs(xDummy3[i][0] - yDummy3[i][0]) < 1e-9);\
    BOOST_CHECK(std::abs(xDummy3[i][1] - yDummy3[i][1]) < 1e-9);\
    BOOST_CHECK(std::abs(xDummy3[i][2] - yDummy3[i][2]) < 1e-9);\
  }

#endif

#endif
