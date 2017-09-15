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
#ifndef __INCLUDED_SINGLESLATER_POP_HPP__
#define __INCLUDED_SINGLESLATER_POP_HPP__

#include <singleslater.hpp>

namespace ChronusQ {

  template <typename T>
  void SingleSlater<T>::populationAnalysis() {

    const size_t NB = aoints.basisSet().nBasis;
    T* SCR  = this->memManager.template malloc<T>(NB*NB);
    T* SCR2 = this->memManager.template malloc<T>(NB*NB);

    // Mulliken population analysis
    mullikenCharges.clear();

    Gemm('N','N',NB,NB,NB,T(1.),aoints.overlap,NB,this->onePDM[SCALAR],NB,
      T(0.),SCR,NB);

    for(auto iAtm = 0; iAtm < aoints.molecule().nAtoms; iAtm++) {

      size_t iEnd;
      if( iAtm == aoints.molecule().nAtoms-1 )
        iEnd = NB;
      else
        iEnd = aoints.basisSet().mapCen2BfSt[iAtm+1];

      size_t iSt = aoints.basisSet().mapCen2BfSt[iAtm];

      mullikenCharges.emplace_back(aoints.molecule().atoms[iAtm].atomicNumber);
      for(auto i = iSt; i < iEnd; i++)
        mullikenCharges.back() -= std::real(SCR[i*(NB+1)]);
    } 


    // Loewdin population analysis
    lowdinCharges.clear();

/*
    Gemm('N','N',NB,NB,NB,T(1.),aoints.ortho1,NB,this->onePDM[SCALAR],NB,
      T(0.),SCR2,NB);
    Gemm('N','C',NB,NB,NB,T(1.),aoints.ortho1,NB,SCR2,NB,T(0.),SCR,NB);
*/

    for(auto iAtm = 0; iAtm < aoints.molecule().nAtoms; iAtm++) {

      size_t iEnd;
      if( iAtm == aoints.molecule().nAtoms-1 )
        iEnd = NB;
      else
        iEnd = aoints.basisSet().mapCen2BfSt[iAtm+1];

      size_t iSt = aoints.basisSet().mapCen2BfSt[iAtm];

      lowdinCharges.emplace_back(aoints.molecule().atoms[iAtm].atomicNumber);
      for(auto i = iSt; i < iEnd; i++)
        lowdinCharges.back() -= std::real(this->onePDMOrtho[SCALAR][i*(NB+1)]);
    } 


    this->memManager.free(SCR,SCR2);


  };

}; // namespace ChronusQ


#endif
