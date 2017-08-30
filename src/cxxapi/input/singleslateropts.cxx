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

    out << "  *** Parsing QM.REFERENCE options ***\n";

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

    std::string RCflag;

    // Determine the Real/Complex flag
    if( tokens.size() == 1 )      RCflag = "AUTO";
    else if( tokens.size() == 2 ) RCflag = tokens[0];
    else CErr("QM.REFERENCE Field not valid",out);


    // Kohn-Sham Keywords
    std::vector<std::string> KSRefs {
      "SLATER", 
      "B88",
      "LSDA",
      "SVWN5",
      "BLYP",
      "PBEXPBEC",
      "B3LYP",
      "PBE0",
      "BHANDHLYP",
      "BHANDH",
    };

    // All refernce keywords
    std::vector<std::string> rawRefs(KSRefs);
    rawRefs.insert(rawRefs.begin(),"HF");

    // Construct R/U/G/X2C reference keywords
    std::vector<std::string> RRefs, URefs, GRefs, X2CRefs;
    for(auto &f : rawRefs) {
      RRefs.emplace_back( "R" + f );
      URefs.emplace_back( "U" + f );
      GRefs.emplace_back( "G" + f );
      X2CRefs.emplace_back( "X2C" + f );
    }


    // This is the reference string to be parsed
    std::string refString = tokens.back();
   

    // Determine type of reference
    bool isRawRef = 
      std::find(rawRefs.begin(),rawRefs.end(),refString) != rawRefs.end();
    bool isRRef   = 
      std::find(RRefs.begin(),RRefs.end(),refString) != RRefs.end();
    bool isURef   = 
      std::find(URefs.begin(),URefs.end(),refString) != URefs.end();
    bool isGRef   = 
      std::find(GRefs.begin(),GRefs.end(),refString) != GRefs.end();
    bool isX2CRef = 
      std::find(X2CRefs.begin(),X2CRefs.end(),refString) != X2CRefs.end();

    // Throw an error if not a valid reference keyword
    if( not isRawRef and not isRRef and not isURef and not isGRef and 
        not isX2CRef )
      CErr(refString + " is not a valid QM.REFERENCE",out);

    // Cleanup the reference string
    if( not isRawRef )
      if( isX2CRef ) refString.erase(0,3);
      else           refString.erase(0,1);


    // Handle KS related queries
    bool isKSRef = 
      std::find(KSRefs.begin(),KSRefs.end(),refString) != KSRefs.end();

    std::string funcName;
    if( isKSRef )
      funcName = refString;

    // Build Functional List
    std::vector<std::shared_ptr<DFTFunctional>> funcList;
    if( isKSRef ) {
      if( not funcName.compare("B88") )
        funcList.push_back(
          std::dynamic_pointer_cast<DFTFunctional>(
            std::make_shared<BEightyEight>()
          )
        );

      if( not funcName.compare("SLATER") )
        funcList.push_back(
          std::dynamic_pointer_cast<DFTFunctional>(
            std::make_shared<SlaterExchange>()
          )
        );

      if( not funcName.compare("LSDA") or not funcName.compare("LDA") ) {

        funcList.push_back(
          std::dynamic_pointer_cast<DFTFunctional>(
            std::make_shared<SlaterExchange>()
          )
        );

        funcList.push_back(
          std::dynamic_pointer_cast<DFTFunctional>(
            std::make_shared<VWNV>()
          )
        );

      }

      if( not funcName.compare("BLYP") ) {

        funcList.push_back(
          std::dynamic_pointer_cast<DFTFunctional>(
            std::make_shared<BEightyEight>()
          )
        );

        funcList.push_back(
          std::dynamic_pointer_cast<DFTFunctional>(
            std::make_shared<LYP>()
          )
        );

      }

      if( not funcName.compare("PBEXPBEC") ) {

        funcList.push_back(
          std::dynamic_pointer_cast<DFTFunctional>(
            std::make_shared<PBEX>()
          )
        );

        funcList.push_back(
          std::dynamic_pointer_cast<DFTFunctional>(
            std::make_shared<PBEC>()
          )
        );

      }

      if( not funcName.compare("B3LYP") ) 
        funcList.push_back(
          std::dynamic_pointer_cast<DFTFunctional>(
            std::make_shared<B3LYP>()
          )
        );

      if( not funcName.compare("PBE0") ) 
        funcList.push_back(
          std::dynamic_pointer_cast<DFTFunctional>(
            std::make_shared<PBE0>()
          )
        );

      if( not funcName.compare("BHANDH") )
        funcList.push_back(
          std::dynamic_pointer_cast<DFTFunctional>(
            std::make_shared<BHANDH>()
          )
        );

      if( not funcName.compare("BHANDHLYP") )
        funcList.push_back(
          std::dynamic_pointer_cast<DFTFunctional>(
            std::make_shared<BHANDHLYP>()
          )
        );

    }
      


    // Setup references
    size_t nC = 1; bool iCS;

    // Raw reference
    if( isRawRef ) {
      out << "  *** Auto-determination of reference: " << refString << " -> ";
      iCS = aoints.molecule().multip == 1;

      if(iCS) out << "R" << refString;
      else    out << "U" << refString;

      out << " ***" << std::endl;
      
    } else if( isRRef )
      if( aoints.molecule().multip != 1 )
        CErr("Spin-Restricted Reference only valid for singlet spin multiplicities",out);
      else
        iCS = true;
    else if( isURef )
      iCS = false;
    else if( isGRef or isX2CRef ) {
      iCS = false; nC = 2;
    }

    // Sanity Checks
    if( nC == 2 and not RCflag.compare("REAL") )
      CErr("Real + Two-Component not valid",out);

    if( nC == 2 and isKSRef )
      CErr("Kohn-Sham + Two-Component not valid",out);

    // Determine Real/Complex if need be
    if(not RCflag.compare("AUTO") ) {
      if( nC == 2 )
        RCflag = "COMPLEX";
      else
        RCflag = "REAL";

      out << "  *** Auto-determination of wave function field: AUTO -> " 
          << RCflag << " ***" << std::endl;
    }

    out << "\n\n";


    // Override core hamiltoninan type for X2C
    if( isX2CRef ) 
      aoints.coreType = EXACT_2C;



    std::shared_ptr<SingleSlaterBase> ss;

    if( not RCflag.compare("REAL") )
      if( isKSRef )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<KohnSham<double>>(
              funcName,funcList,aoints,nC,iCS
            )
          );
      else
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<HartreeFock<double>>(
              aoints,nC,iCS
            )
          );
    else if( not RCflag.compare("COMPLEX") )
      if( isKSRef )
#if 0
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<KohnSham<dcomplex>>(
              aoints,nC,iCS
            )
          );
#else
        CErr();
#endif

      else if( isX2CRef )
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

    return ss;

  }; // CQSingleSlaterOptions

}; // namespace ChronusQ
