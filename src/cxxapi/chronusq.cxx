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

#include <cxxapi/input.hpp>
#include <cxxapi/output.hpp>
#include <cxxapi/options.hpp>
#include <cxxapi/boilerplate.hpp>

#include <util/matout.hpp>
#include <util/threads.hpp>

#include <memmanager.hpp>
#include <cerr.hpp>
#include <molecule.hpp>
#include <basisset.hpp>
#include <aointegrals.hpp>
#include <singleslater.hpp>

#include <cqlinalg/blasext.hpp>

using namespace ChronusQ;

int main(int argc, char *argv[]) {

  ChronusQ::initialize();

  std::string inFileName, outFileName;

  // Parse command line options
  if(argc < 2) { // No Options

    CErr("No Command Line Arguements Given",std::cout);

  } else if(argc == 2) { // Just Input File

    inFileName = argv[1];
    std::vector<std::string> tokens;
    split(tokens,inFileName,".");

    if(tokens.size() == 1) outFileName = inFileName + ".out";
    else outFileName = tokens[0] + ".out";

  } else { // Variable Argc

    int c;
    while((c = getopt(argc,argv,"i:o:")) != -1) {
      switch(c) {
        case('i'):
          inFileName = optarg;
          break;
        case('o'):
          outFileName = optarg;
          break;
        default:
          abort();
      };
    };

  }

  // Create output files
  std::ofstream outfile(outFileName);


  // Redirect output to output file
  std::streambuf *coutbuf = std::cout.rdbuf();
  std::cout.rdbuf(outfile.rdbuf());

  // Output CQ header
  CQOutputHeader(std::cout);

  // Parse Input File
  CQInputFile input(inFileName);

  // Create Molecule and BasisSet objects
  Molecule mol(std::move(CQMoleculeOptions(std::cout,input)));
  BasisSet basis(std::move(CQBasisSetOptions(std::cout,input,mol)));


  CQMemManager memManager(2e9);
  AOIntegrals aoints(memManager,mol,basis);
  auto ss = CQSingleSlaterOptions(std::cout,input,aoints);

  CQSCFOptions(std::cout,input,*ss);
  CQIntsOptions(std::cout,input,aoints);

  aoints.computeCoreHam();
  aoints.computeERI();

  SetNumThreads(1);

  ss->formGuess();
  ss->SCF();

#if 0

  double *SX  = memManager.malloc<double>(basis.nBasis*basis.nBasis);
  double *SX2 = memManager.malloc<double>(basis.nBasis*basis.nBasis);
  double *SX3 = memManager.malloc<double>(basis.nBasis*basis.nBasis);
  double *SX4 = memManager.malloc<double>(basis.nBasis*basis.nBasis);

  double *RX  = memManager.malloc<double>(basis.nBasis*basis.nBasis);
  double *RX2 = memManager.malloc<double>(basis.nBasis*basis.nBasis);
  double *RX3 = memManager.malloc<double>(basis.nBasis*basis.nBasis);
  double *RX4 = memManager.malloc<double>(basis.nBasis*basis.nBasis);

  dcomplex *HX  = memManager.malloc<dcomplex>(basis.nBasis*basis.nBasis);
  dcomplex *HX2 = memManager.malloc<dcomplex>(basis.nBasis*basis.nBasis);
  dcomplex *HX3 = memManager.malloc<dcomplex>(basis.nBasis*basis.nBasis);
  dcomplex *HX4 = memManager.malloc<dcomplex>(basis.nBasis*basis.nBasis);

  dcomplex *ZX  = memManager.malloc<dcomplex>(basis.nBasis*basis.nBasis);
  dcomplex *ZX2 = memManager.malloc<dcomplex>(basis.nBasis*basis.nBasis);
  dcomplex *ZX3 = memManager.malloc<dcomplex>(basis.nBasis*basis.nBasis);
  dcomplex *ZX4 = memManager.malloc<dcomplex>(basis.nBasis*basis.nBasis);

  size_t NB = basis.nBasis;


  std::random_device r;
  std::default_random_engine e(r());
  std::uniform_real_distribution<> dis(-50,68);

  // Create a random real nonsymmetric matrix
  double *Rand = memManager.malloc<double>(NB*NB);

  for(auto i = 0; i < NB; i++)
  for(auto j = 0; j < NB; j++)
    Rand[j + i*NB] = dis(e);



  // Create a random complex nonsymmetric matrix
  dcomplex *ZRand = memManager.malloc<dcomplex>(NB*NB);
  for(auto i = 0; i < NB; i++)
  for(auto j = 0; j < NB; j++)
    ZRand[j + i*NB] = dcomplex(dis(e),dis(e));


  // Symmetrize the complex matrix
  dcomplex *Her = memManager.malloc<dcomplex>(NB*NB);
  std::copy_n(ZRand,NB*NB,Her);
  HerMat('L',NB,Her,NB);


  prettyPrintSmart(std::cout,"Non-Hermetian matrix",ZRand,NB,NB,NB);
  prettyPrintSmart(std::cout,"Hermetian matrix",Her,NB,NB,NB);


  std::vector<TwoBodyContraction<double,double>> cont = 
    { 
      {aoints.overlap, SX , true, COULOMB },
      {aoints.overlap, SX2, true, EXCHANGE} 
    };
  std::vector<TwoBodyContraction<double,double>> cont2 = 
    { 
      {aoints.overlap, SX3, true, COULOMB },
      {aoints.overlap, SX4, true, EXCHANGE} 
    };

  std::vector<TwoBodyContraction<double,double>> cont3 = 
    { 
      {Rand, RX , false, COULOMB },
      {Rand, RX2, false, EXCHANGE} 
    };
  std::vector<TwoBodyContraction<double,double>> cont4 = 
    { 
      {Rand, RX3, false, COULOMB },
      {Rand, RX4, false, EXCHANGE} 
    };




  std::vector<TwoBodyContraction<dcomplex,dcomplex>> cont5 = 
    { 
      {Her, HX , true, COULOMB },
      {Her, HX2, true, EXCHANGE} 
    };
  std::vector<TwoBodyContraction<dcomplex,dcomplex>> cont6 = 
    { 
      {Her, HX3, true, COULOMB },
      {Her, HX4, true, EXCHANGE} 
    };

  std::vector<TwoBodyContraction<dcomplex,dcomplex>> cont7 = 
    { 
      {ZRand, ZX , false, COULOMB },
      {ZRand, ZX2, false, EXCHANGE} 
    };
  std::vector<TwoBodyContraction<dcomplex,dcomplex>> cont8 = 
    { 
      {ZRand, ZX3, false, COULOMB },
      {ZRand, ZX4, false, EXCHANGE} 
    };


  std::fill_n(SX,NB*NB,0.);
  std::fill_n(SX2,NB*NB,0.);
  std::fill_n(SX3,NB*NB,0.);
  std::fill_n(SX4,NB*NB,0.);

  std::fill_n(RX,NB*NB,0.);
  std::fill_n(RX2,NB*NB,0.);
  std::fill_n(RX3,NB*NB,0.);
  std::fill_n(RX4,NB*NB,0.);

  std::fill_n(HX,NB*NB,0.);
  std::fill_n(HX2,NB*NB,0.);
  std::fill_n(HX3,NB*NB,0.);
  std::fill_n(HX4,NB*NB,0.);

  std::fill_n(ZX,NB*NB,0.);
  std::fill_n(ZX2,NB*NB,0.);
  std::fill_n(ZX3,NB*NB,0.);
  std::fill_n(ZX4,NB*NB,0.);


  std::cerr << "Starting Real Symmetric Direct Contraction " << std::endl;
  aoints.twoBodyContractDirect(cont);
//prettyPrintSmart(std::cout,"Real SJ-Direct",SX,NB,NB,NB);
//prettyPrintSmart(std::cout,"Real SK-Direct",SX2,NB,NB,NB);

  std::cerr << "Starting Real Symmetric Incore Contraction " << std::endl;
  aoints.twoBodyContractIncore(cont2);
//prettyPrintSmart(std::cout,"Real SJ-Incore",SX3,NB,NB,NB);
//prettyPrintSmart(std::cout,"Real SK-Incore",SX4,NB,NB,NB);

#if 1
  std::cerr << "Starting Real Nonsymmetric Direct Contraction " << std::endl;
  aoints.twoBodyContractDirect(cont3);
//prettyPrintSmart(std::cout,"Real RJ-Direct",RX,NB,NB,NB);
//prettyPrintSmart(std::cout,"Real RK-Direct",RX2,NB,NB,NB);

  std::cerr << "Starting Real Nonsymmetric Incore Contraction " << std::endl;
  aoints.twoBodyContractIncore(cont4);
//prettyPrintSmart(std::cout,"Real RJ-Incore",RX3,NB,NB,NB);
//prettyPrintSmart(std::cout,"Real RK-Incore",RX4,NB,NB,NB);

  std::cerr << "Starting Complex Hermetian Direct Contraction " << std::endl;
  aoints.twoBodyContractDirect(cont5);
//prettyPrintSmart(std::cout,"Complex HJ-Direct",HX,NB,NB,NB);
  prettyPrintSmart(std::cout,"Complex HK-Direct",HX2,NB,NB,NB);

  std::cerr << "Starting Complex Hermetian Incore Contraction " << std::endl;
  aoints.twoBodyContractIncore(cont6);
//prettyPrintSmart(std::cout,"Complex HJ-Incore",HX3,NB,NB,NB);
  prettyPrintSmart(std::cout,"Complex HK-Incore",HX4,NB,NB,NB);

  std::cerr << "Starting Complex Nonsymmetric Direct Contraction " << std::endl;
  aoints.twoBodyContractDirect(cont7);
//prettyPrintSmart(std::cout,"Complex ZJ-Direct",ZX,NB,NB,NB);
//prettyPrintSmart(std::cout,"Complex ZK-Direct",ZX2,NB,NB,NB);

  std::cerr << "Starting Complex Nonsymmetric Incore Contraction " << std::endl;
  aoints.twoBodyContractIncore(cont8);
//prettyPrintSmart(std::cout,"Complex ZJ-Incore",ZX3,NB,NB,NB);
//prettyPrintSmart(std::cout,"Complex ZK-Incore",ZX4,NB,NB,NB);

#endif



  double maxDiff(0.);
  for(auto i = 0; i < NB*NB; i++)
    maxDiff = std::max(maxDiff,std::abs(SX[i] - SX3[i]));

  std::cerr << "Real Symmetric J = " << maxDiff << std::endl;

  maxDiff = 0.;
  for(auto i = 0; i < NB*NB; i++)
    maxDiff = std::max(maxDiff,std::abs(SX2[i] - SX4[i]));

  std::cerr << "Real Symmetric K = " << maxDiff << std::endl;

  maxDiff =  0.;
  for(auto i = 0; i < NB*NB; i++)
    maxDiff = std::max(maxDiff,std::abs(RX[i] - RX3[i]));

  std::cerr << "Real Nonsymmetric J = " << maxDiff << std::endl;

  maxDiff = 0.;
  for(auto i = 0; i < NB*NB; i++)
    maxDiff = std::max(maxDiff,std::abs(RX2[i] - RX4[i]));

  std::cerr << "Real Nonsymmetric K = " << maxDiff << std::endl;

  maxDiff = 0.;
  for(auto i = 0; i < NB*NB; i++)
    maxDiff = std::max(maxDiff,std::abs(HX[i] - HX3[i]));

  std::cerr << "Complex Hermetian J = " << maxDiff << std::endl;

  maxDiff = 0.;
  for(auto i = 0; i < NB*NB; i++)
    maxDiff = std::max(maxDiff,std::abs(HX2[i] - HX4[i]));

  std::cerr << "Complex Hermetian K = " << maxDiff << std::endl;

  maxDiff =  0.;
  for(auto i = 0; i < NB*NB; i++)
    maxDiff = std::max(maxDiff,std::abs(ZX[i] - ZX3[i]));

  std::cerr << "Complex Nonsymmetric J = " << maxDiff << std::endl;

  maxDiff = 0.;
  for(auto i = 0; i < NB*NB; i++)
    maxDiff = std::max(maxDiff,std::abs(ZX2[i] - ZX4[i]));

  std::cerr << "Complex Nonsymmetric K = " << maxDiff << std::endl;

//prettyPrintSmart(std::cout,"Diff",aoints.overlap,NB,NB,NB);
  memManager.free(SX,SX2,SX3,SX4,Rand,RX,RX2,RX3,RX4,ZRand,HX,HX2,HX3,HX4,Her,ZX,ZX2,ZX3,ZX4);

#endif


  // Output CQ footer
  CQOutputFooter(std::cout);

  // Reset the standard out
  std::cout.rdbuf(coutbuf);

  ChronusQ::finalize();

  return 0;
}

