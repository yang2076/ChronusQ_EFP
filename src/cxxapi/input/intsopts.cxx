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
   *  \brief Optionally set the control parameters for an
   *  AOIntegrals object
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] input  Input file datastructure
   *  \param [in] aoints AOIntegrals object 
   *
   *
   */ 
  void CQIntsOptions(std::ostream &out, CQInputFile &input, AOIntegrals &aoi) {

    // Parse integral algorithm
    std::string ALG = "DIRECT";
    OPTOPT( ALG = input.getData<std::string>("INTS.ALG"); )
    trim(ALG);

    if( not ALG.compare("DIRECT") )
      aoi.cAlg = CONTRACTION_ALGORITHM::DIRECT;
    else if( not ALG.compare("INCORE") )
      aoi.cAlg = CONTRACTION_ALGORITHM::INCORE;
    else
      CErr(ALG + "not a valid INTS.ALG",out);

    
    // Parse Schwartz threshold
    OPTOPT( aoi.threshSchwartz = input.getData<double>("INTS.SCHWARTZ"); )

    out << aoi << std::endl;

  }; // CQIntsOptions


}; // namespace ChronusQ
