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
#include <func.hpp>

#include <cxxapi/input.hpp>
#include <cxxapi/options.hpp>
#include <cxxapi/boilerplate.hpp>

#include <util/threads.hpp>

#include <memmanager.hpp>
#include <cerr.hpp>
#include <molecule.hpp>
#include <basisset.hpp>
#include <aointegrals.hpp>

#include <cqlinalg/blasext.hpp>


using namespace ChronusQ;

#ifdef _CQ_GENERATE_TESTS

   /**
    *  \brief Global Boost.Test fixture for Contraction
    *  tests in the case of test generation.
    *
    *  Creates / overwrites the reference file for the contraction
    *  data.
    *
    */ 
  struct ContractConfig {

    ContractConfig() { 
      SafeFile refFile(FUNC_REFERENCE "contract.hdf5"); 
      refFile.createFile();
    }

  };

  BOOST_GLOBAL_FIXTURE( ContractConfig );

#endif





/**
 *  \brief Generate a randon number of different
 *  types using the STL random number generator
 */
template <typename T>
T RAND_NUMBER(std::default_random_engine &e, 
  std::uniform_real_distribution<> &dis);

template <>
double RAND_NUMBER(std::default_random_engine &e, 
  std::uniform_real_distribution<> &dis) {

  return dis(e);

}; // RAND_NUMBER<double>

template <>
dcomplex RAND_NUMBER(std::default_random_engine &e, 
  std::uniform_real_distribution<> &dis) {

  return dcomplex(dis(e),dis(e));

}; // RAND_NUMBER<dcomplex>


/**
 *  Structure of matricies for contration tests
 */
enum HER_MAT {
  HERMETIAN,
  NONHERMETIAN
};

/**
 *  \brief Template wrapper around functions to optionally
 *  symmetrize a matrix. Null call for H == NONHERMETIAN
 */
template <typename T, HER_MAT H>
void HerMatLoc(char UPLO, size_t N, T *A) {

  if(H == HERMETIAN) ChronusQ::HerMat(UPLO,N,A,N);
  else { ; }

};

// Set up the integrals/memory for a direct contraction test
#define CONTRACT_BUILD(FIELD,HER,TYPE) \
  CQInputFile input(FUNC_INPUT "contract_ref.inp");\
  \
  auto memManager = CQMiscOptions(std::cout,input); \
  \
  Molecule mol(std::move(CQMoleculeOptions(std::cout,input))); \
  BasisSet basis(std::move(CQBasisSetOptions(std::cout,input,mol))); \
  AOIntegrals aoints(*memManager,mol,basis); \
  \
  size_t NB = basis.nBasis; \
  FIELD *SX  = memManager->malloc<FIELD>(NB*NB); \
  FIELD *SX2 = memManager->malloc<FIELD>(NB*NB); \
  FIELD *Rand = memManager->malloc<FIELD>(NB*NB); \
  std::vector<TwoBodyContraction<FIELD,FIELD>> cont = \
    { { Rand, SX, (HER == HERMETIAN), TYPE  } };\
  std::fill_n(SX,NB*NB,0.);\
  \





#ifdef _CQ_GENERATE_TESTS


// If we're generating the tests, generate the random matricies
// and perform the contraction incore
#define CONTRACT_TEST(FIELD,HER,TYPE,STORAGE) \
  /* Get the reference File */ \
  SafeFile refFile(FUNC_REFERENCE "contract.hdf5"); \
  \
  CONTRACT_BUILD(FIELD,HER,TYPE) \
  \
  aoints.computeERI(); \
  \
  std::random_device r; \
  std::default_random_engine e(r());\
  std::uniform_real_distribution<> dis(-50,68); \
  \
  for(auto i = 0; i < NB; i++)\
  for(auto j = 0; j < NB; j++)\
    Rand[j + i*NB] = RAND_NUMBER<FIELD>(e,dis);\
  \
  HerMatLoc<FIELD,HER>('U',NB,Rand); \
  refFile.safeWriteData(STORAGE "/X",Rand,{NB,NB});\
  \
  aoints.twoBodyContractIncore(cont);\
  refFile.safeWriteData(STORAGE "/AX",SX,{NB,NB});\
  memManager->free(SX,SX2,Rand);


#else


// If we're not generating the tests, read the reference data off
// disk and perform / compare the contraction incore
#define CONTRACT_TEST(FIELD,HER,TYPE,STORAGE) \
  SafeFile refFile(FUNC_REFERENCE "contract.hdf5",true);\
  \
  CONTRACT_BUILD(FIELD,HER,TYPE) \
  \
  refFile.readData(STORAGE "/X",Rand);\
  refFile.readData(STORAGE "/AX",SX2);\
  \
  aoints.twoBodyContractDirect(cont);\
  \
  double maxDiff(0.);\
  for(auto i = 0; i < NB*NB; i++) \
    maxDiff = std::max(maxDiff,std::abs(SX[i] - SX2[i]));\
  \
  BOOST_CHECK(maxDiff < 1e-10);\
  memManager->free(SX,SX2,Rand);



#endif


// Direct contract test suite
BOOST_AUTO_TEST_SUITE( DIRECT_CONTRACTION )



// Real contraction test suite
BOOST_AUTO_TEST_SUITE( REAL_DIRECT_CONTRACTION )

// Real Hermetian "J" contraction test
BOOST_FIXTURE_TEST_CASE( HER_J_CONTRACT, SerialJob ) {

  CONTRACT_TEST(double,HERMETIAN,COULOMB,"CONTRACTION/HER/REAL/J");

}

// Real Non-Hermetian "J" contraction test
BOOST_FIXTURE_TEST_CASE( NONHER_J_CONTRACT, SerialJob ) {

  CONTRACT_TEST(double,NONHERMETIAN,COULOMB,"CONTRACTION/NONHER/REAL/J");

}

// Real Hermetian "K" contraction test
BOOST_FIXTURE_TEST_CASE( HER_K_CONTRACT, SerialJob ) {

  CONTRACT_TEST(double,HERMETIAN,EXCHANGE,"CONTRACTION/HER/REAL/K");

}

// Real Non-Hermetian "K" contraction test
BOOST_FIXTURE_TEST_CASE( NONHER_K_CONTRACT, SerialJob ) {

  CONTRACT_TEST(double,NONHERMETIAN,COULOMB,"CONTRACTION/NONHER/REAL/K");

}

// Parallel Real direct contraction tests
#ifdef _CQ_DO_PARTESTS

// Parallel Real Hermetian "J" contraction test
BOOST_FIXTURE_TEST_CASE( PAR_HER_J_CONTRACT, ParallelJob ) {

  CONTRACT_TEST(double,HERMETIAN,COULOMB,"CONTRACTION/HER/REAL/J");

}

// Parallel Real Non-Hermetian "J" contraction test
BOOST_FIXTURE_TEST_CASE( PAR_NONHER_J_CONTRACT, ParallelJob ) {

  CONTRACT_TEST(double,NONHERMETIAN,COULOMB,"CONTRACTION/NONHER/REAL/J");

}

// Parallel Real Hermetian "K" contraction test
BOOST_FIXTURE_TEST_CASE( PAR_HER_K_CONTRACT, ParallelJob ) {

  CONTRACT_TEST(double,HERMETIAN,EXCHANGE,"CONTRACTION/HER/REAL/K");

}

// Parallel Real Non-Hermetian "K" contraction test
BOOST_FIXTURE_TEST_CASE( PAR_NONHER_K_CONTRACT, ParallelJob ) {

  CONTRACT_TEST(double,NONHERMETIAN,COULOMB,"CONTRACTION/NONHER/REAL/K");

}


#endif

// End real contraction test suite
BOOST_AUTO_TEST_SUITE_END()


// Complex contraction test suite
BOOST_AUTO_TEST_SUITE( COMPLEX_DIRECT_CONTRACTION )

// Complex Hermetian "J" contraction test
BOOST_FIXTURE_TEST_CASE( HER_J_CONTRACT, SerialJob ) {

  CONTRACT_TEST(dcomplex,HERMETIAN,COULOMB,"CONTRACTION/HER/COMPLEX/J");

}

// Complex Non-Hermetian "J" contraction test
BOOST_FIXTURE_TEST_CASE( NONHER_J_CONTRACT, SerialJob ) {

  CONTRACT_TEST(dcomplex,NONHERMETIAN,COULOMB,"CONTRACTION/NONHER/COMPLEX/J");

}


// Complex Hermetian "K" contraction test
BOOST_FIXTURE_TEST_CASE( HER_K_CONTRACT, SerialJob ) {

  CONTRACT_TEST(dcomplex,HERMETIAN,EXCHANGE,"CONTRACTION/HER/COMPLEX/K");

}

// Complex Non-Hermetian "K" contraction test
BOOST_FIXTURE_TEST_CASE( NONHER_K_CONTRACT, SerialJob ) {

  CONTRACT_TEST(dcomplex,HERMETIAN,EXCHANGE,"CONTRACTION/NONHER/COMPLEX/K");

}



// Parallel Complex direct contraction tests
#ifdef _CQ_DO_PARTESTS 

// Parallel Complex Hermetian "J" contraction test
BOOST_FIXTURE_TEST_CASE( PAR_HER_J_CONTRACT, ParallelJob ) {

  CONTRACT_TEST(dcomplex,HERMETIAN,COULOMB,"CONTRACTION/HER/COMPLEX/J");

}

// Parallel Complex Non-Hermetian "J" contraction test
BOOST_FIXTURE_TEST_CASE( PAR_NONHER_J_CONTRACT, ParallelJob ) {

  CONTRACT_TEST(dcomplex,NONHERMETIAN,COULOMB,"CONTRACTION/NONHER/COMPLEX/J");

}


// Parallel Complex Hermetian "K" contraction test
BOOST_FIXTURE_TEST_CASE( PAR_HER_K_CONTRACT, ParallelJob ) {

  CONTRACT_TEST(dcomplex,HERMETIAN,EXCHANGE,"CONTRACTION/HER/COMPLEX/K");

}

// Parallel Complex Non-Hermetian "K" contraction test
BOOST_FIXTURE_TEST_CASE( PAR_NONHER_K_CONTRACT, ParallelJob ) {

  CONTRACT_TEST(dcomplex,HERMETIAN,EXCHANGE,"CONTRACTION/NONHER/COMPLEX/K");

}

#endif


// End complex contraction suite
BOOST_AUTO_TEST_SUITE_END()

// End direct contraction suite
BOOST_AUTO_TEST_SUITE_END()
