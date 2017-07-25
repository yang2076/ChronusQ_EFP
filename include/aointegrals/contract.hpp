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
#ifndef __INCLUDED_AOINTEGRALS_CONTRACT_HPP__
#define __INCLUDED_AOINTEGRALS_CONTRACT_HPP__


#include <aointegrals.hpp>
#include <cqlinalg/blas3.hpp>

// Use stupid but bullet proof incore contraction for debug
//#define _BULLET_PROOF_INCORE

namespace ChronusQ {

  /**
   *  \brief Perform various tensor contractions of the full ERI
   *  tensor in core.
   *
   *  Currently supports
   *    - Coulomb-type (34,12) contractions
   *    - Exchange-type (23,12) contractions
   *
   *  Works with both real and complex matricies
   *
   *  \param [in/out] list Contains information pertinent to the 
   *    matricies to be contracted with. See TwoBodyContraction
   *    for details
   */ 
  template <typename T>
  void AOIntegrals::twoBodyContractIncore(
    std::vector<TwoBodyContraction<T>> &list) {

    size_t NB3 = basisSet_.nBasis * nSQ_;

    // Loop over matricies to contract with
    for(auto &C : list) {
      // Sanity check of dimensions
      assert(memManager_.template getSize<T>(C.X)  == nSQ_);
      assert(memManager_.template getSize<T>(C.AX) == nSQ_);

      // Coulomb-type (34,12) ERI contraction
      // AX(mn) = (mn | kl) X(kl)
      if( C.contType == COULOMB ) {

        #ifdef _BULLET_PROOF_INCORE

        for(auto i = 0; i < basisSet_.nBasis; ++i)
        for(auto j = 0; j < basisSet_.nBasis; ++j)
        for(auto k = 0; k < basisSet_.nBasis; ++k)
        for(auto l = 0; l < basisSet_.nBasis; ++l) {
          C.AX[i + j*basisSet_.nBasis] +=
            ERI[i + j*basisSet_.nBasis + k*nSQ_ + l*NB3] *
            C.X[k + l*basisSet_.nBasis];
        }

        #else
  
        Gemm('N','N',nSQ_,1,nSQ_,T(1.),ERI,nSQ_,C.X,nSQ_,T(0.),C.AX,nSQ_);

        #endif


      // Exchange-type (23,12) ERI contraction
      // AX(mn) = (mk |ln) X(kl)
      } else if( C.contType == EXCHANGE ) {

        #ifdef _BULLET_PROOF_INCORE

        for(auto i = 0; i < basisSet_.nBasis; ++i)
        for(auto j = 0; j < basisSet_.nBasis; ++j)
        for(auto k = 0; k < basisSet_.nBasis; ++k)
        for(auto l = 0; l < basisSet_.nBasis; ++l) {
          C.AX[i + j*basisSet_.nBasis] +=
            ERI[i + k*basisSet_.nBasis + l*nSQ_ + j*NB3] *
            C.X[k + l*basisSet_.nBasis];
        }

        #else

        for(auto nu = 0; nu < basisSet_.nBasis; nu++) {
          Gemm('N','N',basisSet_.nBasis,1,nSQ_,T(1.),
            ERI  + nu * NB3,              basisSet_.nBasis,
            C.X,                          nSQ_,T(0.),
            C.AX + nu * basisSet_.nBasis, basisSet_.nBasis);
        }

        #endif
      }
    } // loop over matricies
    
  }; // AOIntegrals::twoBodyContractIncore

}; // namespace ChronusQ

#endif
