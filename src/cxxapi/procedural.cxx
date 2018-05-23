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
#include <cxxapi/procedural.hpp>

#include <util/matout.hpp>
#include <util/threads.hpp>

#include <memmanager.hpp>
#include <cerr.hpp>
#include <molecule.hpp>
#include <basisset.hpp>
#include <aointegrals.hpp>
#include <singleslater.hpp>
#include <realtime.hpp>

#include <cqlinalg/blasext.hpp>

#include <util/files.hpp>

namespace ChronusQ {

  void RunChronusQ(std::string inFileName,
    std::string outFileName, std::string rstFileName,
    std::string scrFileName) {

    // Create the restart and scratch files
    SafeFile rstFile(rstFileName);
    //SafeFile scrFile(scrFileName);

    rstFile.createFile();


    // Redirect output to output file if not STDOUT
    std::shared_ptr<std::ofstream> outfile;
    std::streambuf *coutbuf = std::cout.rdbuf();

    if( outFileName.compare("STDOUT") ) {

      outfile = std::make_shared<std::ofstream>(outFileName);
      std::cout.rdbuf(outfile->rdbuf());

    }


    // Output CQ header
    CQOutputHeader(std::cout);

    // Parse Input File
    CQInputFile input(inFileName);



    // Dump contents of input file into output file
    std::cout << "\n\n\n";
    std::cout << "Input File:\n" << BannerTop << std::endl;
    std::ifstream inStream(inFileName);
    std::istreambuf_iterator<char> begin_src(inStream);
    std::istreambuf_iterator<char> end_src;
    std::ostreambuf_iterator<char> begin_dest(std::cout);
    std::copy(begin_src,end_src,begin_dest);
    inStream.close();
    std::cout << BannerEnd << "\n\n\n" << std::endl;




    // Determine JOB type
    std::string jobType;
    
    try {
      jobType = input.getData<std::string>("QM.JOB");
    } catch (...) {
      CErr("Must Specify QM.JOB",std::cout);
    }


    auto memManager = CQMiscOptions(std::cout,input);


    // Create Molecule and BasisSet objects
    Molecule mol(std::move(CQMoleculeOptions(std::cout,input)));
    BasisSet basis(std::move(CQBasisSetOptions(std::cout,input,mol)));


    AOIntegrals aoints(*memManager,mol,basis);
    auto ss = CQSingleSlaterOptions(std::cout,input,aoints);


    ss->savFile    = rstFile;
    aoints.savFile = rstFile;

    // EM Perturbation for SCF
    EMPerturbation SCFpert;

    CQSCFOptions(std::cout,input,*ss,SCFpert);
    CQIntsOptions(std::cout,input,aoints);


    if( not jobType.compare("SCF") or not jobType.compare("RT") ) {

      aoints.computeCoreHam();

      // If INCORE, compute and store the ERIs
      if(aoints.cAlg == INCORE) aoints.computeERI();

      ss->formGuess();
      ss->SCF(SCFpert);
    }

    if( not jobType.compare("RT") ) {
      auto rt = CQRealTimeOptions(std::cout,input,ss);
      rt->savFile = rstFile;
      rt->doPropagation();
    }

    // Output CQ footer
    CQOutputFooter(std::cout);

    // Reset std::cout
    if(outfile) std::cout.rdbuf(coutbuf);

  }; // RunChronusQ


 }; // namespace ChronusQ
