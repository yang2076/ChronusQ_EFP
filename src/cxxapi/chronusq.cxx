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


#include <cxxapi/procedural.hpp>
#include <cxxapi/boilerplate.hpp>

#include <cerr.hpp>

using namespace ChronusQ;

int main(int argc, char *argv[]) {

  ChronusQ::initialize();

  std::string inFileName, outFileName;
  std::string rstFileName, scrFileName;

  // Parse command line options
  if(argc < 2) { // No Options

    CErr("No Command Line Arguements Given",std::cout);

  } else if(argc == 2) { // Just Input File

    inFileName = argv[1];
    std::vector<std::string> tokens;
    split(tokens,inFileName,".");

    outFileName = tokens[0] + ".out";
    rstFileName = tokens[0] + ".bin";

  // FIXME: Need to fix this for generality in specification
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



  RunChronusQ(inFileName,outFileName,rstFileName,scrFileName);

  ChronusQ::finalize();

  return 0;
}

