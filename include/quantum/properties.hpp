/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2018 Li Research Group (University of Washington)
 *  
 *  This program is free software; you ca redistribute it and/or modify
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
#ifndef __INCLUDED_QUANTUM_PROPERTIES_HPP__
#define __INCLUDED_QUANTUM_PROPERTIES_HPP__

#include <quantum.hpp>

namespace ChronusQ {

  /**
   *  \brief Perform matrix trace with proper component of the density for
   *  property evaluation. Return zeros when necessacy.
   */ 
  template <typename MatsT>
  template <typename RetTyp, DENSITY_TYPE DenTyp, typename Op>
  RetTyp Quantum<MatsT>::OperatorSpinCombine(const Op &op) {
    double rZero(0.);
    bool isReal = std::is_same<double,MatsT>::value;

    size_t DSize = memManager.template getSize(onePDM[SCALAR]);

    // Scalar trace option always valid

    // Zero eval if 1C, closed shell and not scalar trace
    bool isZero = (nC == 1 and iCS) and (DenTyp != DENSITY_TYPE::SCALAR);

    // Zero eval if 1C, open shell and MX / MY trace
    isZero = isZero or ( 
      (nC == 1 and not iCS) and
      (DenTyp == DENSITY_TYPE::MY or DenTyp == DENSITY_TYPE::MX) 
    );

    // Zero eval is 2C, real and MY trace
    isZero = isZero or (
      (nC > 1) and (DenTyp == DENSITY_TYPE::MY and isReal)
    );

    // Catch zero evals
    if(isZero)
      return reinterpret_cast<RetTyp(&)[2]>(rZero)[0];

    // Sanity checks
    assert( memManager.template getSize(op) == DSize );

    // Perform proper trace
    return OperatorTrace<RetTyp>(DSize,onePDM[DenTyp],op);

  }; // Quantum<T>::OperatorSpinCombine


}; // namespace ChronusQ

#endif
