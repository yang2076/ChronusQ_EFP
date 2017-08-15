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
#include <cxxapi/options.hpp>
#include <cerr.hpp>

namespace ChronusQ {

  /**
   *  \brief Construct a SingleSlater object using the input 
   *  file.
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] input  Input file datastructure
   *  \param [in] aoints AOIntegrals object for SingleSlater
   *                     construction
   *
   *  \returns shared_ptr to a SingleSlaterBase object
   *    constructed from the input options.
   *
   */ 
  std::shared_ptr<SingleSlaterBase> CQSingleSlaterOptions(
    std::ostream &out, CQInputFile &input, 
    AOIntegrals &aoints) {


    // Attempt to find reference
    std::string reference;
    try { 
      reference = input.getData<std::string>("QM.REFERENCE");
    } catch(...) {
      CErr("QM.REFERENCE Keyword not found!",out);
    }

    // Digest reference string
    // Trim Spaces
    trim(reference);

    // Split into tokens
    std::vector<std::string> tokens;
    split(tokens,reference);
    for(auto &X : tokens) trim(X);

    std::string RCflag, refString;

    // Determine the Real/Complex flag
    if( tokens.size() == 1 )      RCflag = "AUTO";
    else if( tokens.size() == 2 ) RCflag = tokens[0];
    else CErr("QM.REFERENCE Field not valid",out);

    refString = tokens.back();

    size_t nC = 1; bool iCS;
    if( not refString.compare("HF") ) {
      out << "  *** Auto-determination of reference: HF -> ";
      iCS = aoints.molecule().multip == 1;

      if(iCS) out << "RHF";
      else    out << "UHF";

      out << " ***" << std::endl;
      
    } else if( not refString.compare("RHF") )
      if( aoints.molecule().multip != 1 )
        CErr("Spin-Restricted Reference only valid for singlet spin multiplicities",out);
      else
        iCS = true;
    else if( not refString.compare("UHF") )
      iCS = false;
    else if( not refString.compare("GHF") or not refString.compare("X2C") ) {
      iCS = false; nC = 2;
    }

    if( nC == 2 and not RCflag.compare("REAL") )
      CErr("Real + Two-Component not valid",out);

    // Determine Real/Complex if need be
    if(not RCflag.compare("AUTO") ) {
      if( nC == 2 )
        RCflag = "COMPLEX";
      else
        RCflag = "REAL";

      out << "  *** Auto-determination of wave function field: AUTO -> " 
          << RCflag << " ***" << std::endl;
    }

    // Override core hamiltoninan type for X2C
    if( not refString.compare("X2C") ) 
      aoints.coreType = EXACT_2C;


     bool isHF  = true;
     bool isX2C = not refString.compare("X2C");


    std::shared_ptr<SingleSlaterBase> ss;

    if( isHF ) {
      if( not RCflag.compare("REAL") )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<HartreeFock<double>>(
              aoints,nC,iCS
            )
          );
      else if( not RCflag.compare("COMPLEX") )
        if( isX2C )
          ss = std::dynamic_pointer_cast<SingleSlaterBase>(
              std::make_shared<HartreeFock<dcomplex>>(
                "Exact Two Component","X2C",aoints,nC,iCS
              )
            );
        else
          ss = std::dynamic_pointer_cast<SingleSlaterBase>(
              std::make_shared<HartreeFock<dcomplex>>(
                aoints,nC,iCS
              )
            );
    }


    return ss;

  }; // CQSingleSlaterOptions

}; // namespace ChronusQ
