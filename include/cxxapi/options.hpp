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
#ifndef __INCLUDED_OPTIONS_HPP__
#define __INCLUDED_OPTIONS_HPP__

#include <cxxapi/input.hpp>
#include <cxxapi/procedural.hpp>

#include <molecule.hpp>
#include <basisset.hpp>
#include <aointegrals.hpp>
#include <singleslater.hpp>
#include <realtime.hpp>
#include <response.hpp>
#include <memmanager.hpp>
#include <chronusqefp.hpp>

// Preprocessor directive to aid the digestion of optional 
// input arguments
#define OPTOPT(x) try{ x; } catch(...) { ; }

namespace ChronusQ {

  // Function definitions ofr option parsing. 
  // See src/cxxapi/input/*opts.cxx for documentation

  // Parse the options relating to the Molecule object
  Molecule CQMoleculeOptions(std::ostream &, CQInputFile &);

  void CQMOLECULE_VALID(std::ostream&, CQInputFile &);
  // Parse the options relating to the BasisSet
  BasisSet CQBasisSetOptions(std::ostream &, CQInputFile &, 
    Molecule &);

  void CQBASIS_VALID(std::ostream&, CQInputFile &);

  // Parse the options relating to the SingleSlater 
  // (and variants)
  std::shared_ptr<SingleSlaterBase> CQSingleSlaterOptions(
    std::ostream &, CQInputFile &, std::shared_ptr<AOIntegralsBase> );

  void CQQM_VALID(std::ostream&, CQInputFile &);
  void CQDFTINT_VALID(std::ostream&, CQInputFile &);

  // Parse RT options
  std::shared_ptr<RealTimeBase> CQRealTimeOptions(
    std::ostream &, CQInputFile &, std::shared_ptr<SingleSlaterBase> &
  );

  void CQRT_VALID(std::ostream&, CQInputFile &);

  // Parse Response options
  std::shared_ptr<ResponseBase> CQResponseOptions(
    std::ostream &, CQInputFile &, std::shared_ptr<SingleSlaterBase>
  );

  void CQRESPONSE_VALID(std::ostream&, CQInputFile &);
  void CQMOR_VALID(std::ostream&, CQInputFile &);

  // Parse integral options
  std::shared_ptr<AOIntegralsBase> CQIntsOptions(std::ostream &, 
    CQInputFile &, CQMemManager &, Molecule &, BasisSet &); 

  void CQINTS_VALID(std::ostream&, CQInputFile &);

  // Parse the SCF options
  void CQSCFOptions(std::ostream&, CQInputFile&,
    SingleSlaterBase &, EMPerturbation &);

  void CQSCF_VALID(std::ostream&, CQInputFile &);

  std::shared_ptr<CQMemManager> CQMiscOptions(std::ostream &,
    CQInputFile &);

  void CQMISC_VALID(std::ostream&, CQInputFile &);

  // Parse the EFP options
  void CQEFP_VALID(std::ostream&, CQInputFile&);

  ::efp_opts CQEFPOptions(std::ostream&, CQInputFile&);

  void CQFragCoord(std::ostream&, CQInputFile&, CQMemManager&, struct Fragment_ifm*);


  std::string CQEFPFilePath(std::ostream&, CQInputFile&);  
  std::shared_ptr<EFPBase> CQEFPControl(std::ostream&, CQInputFile&, std::shared_ptr<AOIntegralsBase>&, std::shared_ptr<SingleSlaterBase>&,bool);

  inline void CQINPUT_VALID(std::ostream &out, CQInputFile &input) {

    CQMOLECULE_VALID(out,input);
    CQBASIS_VALID(out,input);
    CQINTS_VALID(out,input);
    CQEFP_VALID(out,input);
    CQQM_VALID(out,input);
    CQDFTINT_VALID(out,input);
    CQSCF_VALID(out,input);
    CQRT_VALID(out,input);
    CQRESPONSE_VALID(out,input);
    CQMOR_VALID(out,input);
    CQMISC_VALID(out,input);

  }


};




#endif
