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

#include <cxxapi/input.hpp>
#include <cxxapi/output.hpp>
#include <cxxapi/options.hpp>
#include <cxxapi/boilerplate.hpp>
#include <cxxapi/procedural.hpp>

#include <util/matout.hpp>
#include <util/threads.hpp>
#include <util/mpi.hpp>
#include <util/files.hpp>

#include <memmanager.hpp>
#include <cerr.hpp>
#include <molecule.hpp>
#include <basisset.hpp>
#include <aointegrals.hpp>
#include <singleslater.hpp>
#include <response.hpp>
#include <realtime.hpp>
#include <itersolver.hpp>
#include <chronusqefp.hpp>

#include <cqlinalg/blasext.hpp>
#include <cqlinalg/eig.hpp>

#include <physcon.hpp>


//#include <cubegen.hpp>

namespace ChronusQ {
#ifdef ENABLE_BCAST_COUNTER
  int bcastCounter = 0;
#endif
  void RunChronusQ(std::string inFileName,
    std::string outFileName, std::string rstFileName,
    std::string scrFileName) {
    int rank = MPIRank();
    int size = MPISize();
    
    // Redirect output to output file if not STDOUT
    std::shared_ptr<std::ofstream> outfile;
    std::streambuf *coutbuf = std::cout.rdbuf();

    if( outFileName.compare("STDOUT") and (rank == 0) ) {
      outfile = std::make_shared<std::ofstream>(outFileName);
      std::cout.rdbuf(outfile->rdbuf());

    }


    // Setup MPI rank files

    std::shared_ptr<std::ofstream> rankfile;
    std::streambuf *cerrbuf = std::cerr.rdbuf();

    if( size > 1 ) {
      
      std::string rankFileName = outFileName + ".mpi." + 
        std::to_string(rank);

      rankfile = std::make_shared<std::ofstream>(rankFileName);
      std::cerr.rdbuf(rankfile->rdbuf());

      std::cerr << "Hello from RANK = " << rank << " / SIZE = " << size 
                << "\n\n";

    }


    std::ostream &output = (rank == 0) ? std::cout : std::cerr;


    // Output CQ header
    CQOutputHeader(output);
    if(rankfile and rank == 0) CQOutputHeader(std::cerr);


    // Parse Input File
    CQInputFile input(inFileName);


    CQINPUT_VALID(output,input);

    // Dump contents of input file into output file
    if( rank == 0 ) {
      std::cout << "\n\n\n";
      std::cout << "Input File:\n" << BannerTop << std::endl;
      std::ifstream inStream(inFileName);
      std::istreambuf_iterator<char> begin_src(inStream);
      std::istreambuf_iterator<char> end_src;
      std::ostreambuf_iterator<char> begin_dest(std::cout);
      std::copy(begin_src,end_src,begin_dest);
      inStream.close();
      std::cout << BannerEnd << "\n\n\n" << std::endl;
    }




    // Determine JOB type
    std::string jobType;
    
    try {
      jobType = input.getData<std::string>("QM.JOB");
    } catch (...) {
      CErr("Must Specify QM.JOB",output);
    }


    auto memManager = CQMiscOptions(output,input);


    // Create Molecule and BasisSet objects
    Molecule mol(std::move(CQMoleculeOptions(output,input)));
    BasisSet basis(std::move(CQBasisSetOptions(output,input,mol)));


    auto aoints = CQIntsOptions(output,input,*memManager,mol,basis);
    auto ss = CQSingleSlaterOptions(output,input,aoints);

    // EM Perturbation for SCF
    EMPerturbation emPert;
    // SCF options
    CQSCFOptions(output,input,*ss,emPert);

    // EFP option for SCF
    bool EFP_bool = (ss->scfControls).EFP_bool;
    auto EFP_ = CQEFPControl(output,input,aoints,ss,EFP_bool);
    auto EFP_1 = EFP_.get();

   

    bool rstExists = false;
    if( ss->scfControls.guess == READMO or 
        ss->scfControls.guess == READDEN ) 
      rstExists = true;

    // Create the restart and scratch files
    SafeFile rstFile(rstFileName, rstExists);
    //SafeFile scrFile(scrFileName);

    if( not rstExists and rank == 0 ) rstFile.createFile();


    if( rank == 0 ) {
      ss->savFile     = rstFile;
      aoints->savFile = rstFile;
    }
      

    if( not jobType.compare("SCF") or not jobType.compare("RT") or 
        not jobType.compare("RESP") ) {

      ss->formCoreH(emPert);

      // If INCORE, compute and store the ERIs
      if(aoints->contrAlg == INCORE) aoints->computeERI(emPert);

      ss->formGuess(EFP_1,EFP_bool);
      ss->SCF(emPert,EFP_1,EFP_bool);
    }


    if( not jobType.compare("RT") ) {

      if( MPISize() > 1 ) CErr("RT + MPI NYI!",output);

      auto rt = CQRealTimeOptions(output,input,ss);
      rt->savFile = rstFile;
      rt->EFP_user_data_change(EFP_1,EFP_bool);
      rt->doPropagation(EFP_1,EFP_bool);

    }

    if( not jobType.compare("RESP") ) {
      std::cout << "Before response" << std::endl;
      std::cout << memManager << std::endl;

      auto resp = CQResponseOptions(output,input,ss);
      resp->savFile = rstFile;
      resp->run();

      if( MPIRank(MPI_COMM_WORLD) == 0 ) resp->printResults(output);
      MPI_Barrier(MPI_COMM_WORLD);

    }

    // Output CQ footer
    CQOutputFooter(output);

    // Reset std::cout and std::cerr
    if(outfile)  std::cout.rdbuf(coutbuf);
    if(rankfile) std::cerr.rdbuf(cerrbuf);

  }; // RunChronusQ


 }; // namespace ChronusQ
