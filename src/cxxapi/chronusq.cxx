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

#include <memmanager.hpp>
#include <cerr.hpp>
#include <molecule.hpp>
#include <basisset.hpp>
#include <aointegrals.hpp>
#include <singleslater.hpp>

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


  CQMemManager memManager(100e6);
  AOIntegrals aoints(memManager,mol,basis);

  aoints.computeAOOneE();
  aoints.computeERI();

/*
  double *SX  = memManager.malloc<double>(basis.nBasis*basis.nBasis);
  double *SX2 = memManager.malloc<double>(basis.nBasis*basis.nBasis);
  std::vector<TwoBodyContraction<double>> cont = 
    { 
      {aoints.overlap, SX , true, COULOMB },
      {aoints.overlap, SX2, true, EXCHANGE} 
    };

  aoints.twoBodyContract(cont);

  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> SMap(aoints.overlap,basis.nBasis,basis.nBasis);
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> SXMap(SX,basis.nBasis,basis.nBasis);
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> SX2Map(SX2,basis.nBasis,basis.nBasis);

  std::cout << SMap << std::endl << std::endl;;
  std::cout << SXMap << std::endl << std::endl;;
  std::cout << SX2Map << std::endl;

  memManager.free(SX,SX2);
*/

  SingleSlater<dcomplex> ss(aoints,1);
  ss.formGuess();
  ss.SCF();


  // Output CQ footer
  CQOutputFooter(std::cout);

  ChronusQ::finalize();

  return 0;
}

