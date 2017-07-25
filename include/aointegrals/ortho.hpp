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
#ifndef __INCLUDED_AOINTEGRALS_ORTHO_HPP__
#define __INCLUDED_AOINTEGRALS_ORTHO_HPP__

#include <aointegrals.hpp>
#include <cqlinalg/blas3.hpp>


namespace ChronusQ {

  /**
   *  \brief Performs the transformation \f$ A' = O_1 A O_1^T\f$ for
   *  a general matrix \f$A\f$ 
   *
   *  \f$ O_1 \f$ is the orthonormalization matrix stored in 
   *  AOIntegrals::ortho1 (see AOIntegrals::ORTHO_TYPE and
   *  AOIntegrals::computeOrtho for details).
   *
   */ 
  template <typename T> 
  void AOIntegrals::Ortho1Trans(T* A, T* TransA) {

    // Make sure that the incoming matricies are of the right size
    assert(memManager_.template getSize(A) == nSQ_);
    assert(memManager_.template getSize(TransA) == nSQ_);

    // Allocate scratch space
    T* SCR = memManager_.template malloc<T>(nSQ_);

    // If necessary, create a complex copy of ortho1
    T* O1 = reinterpret_cast<T*>(ortho1);
    if(std::is_same<T,dcomplex>::value) {
      O1 = memManager_.template malloc<T>(nSQ_);
      std::copy_n(ortho1,nSQ_,O1);
    }

    // Perform transformation
    Gemm('N', 'N', basisSet_.nBasis, basisSet_.nBasis, basisSet_.nBasis, T(1.),
      O1, basisSet_.nBasis, A, basisSet_.nBasis, T(0.), SCR, basisSet_.nBasis);
    Gemm('N', 'T', basisSet_.nBasis, basisSet_.nBasis, basisSet_.nBasis, T(1.),
      SCR, basisSet_.nBasis, O1, basisSet_.nBasis, T(0.), TransA,
      basisSet_.nBasis);

    // Free up scratch space
    memManager_.free(SCR);
    if(std::is_same<T,dcomplex>::value) memManager_.free(O1);
    
  }; // AOIntegrals::Ortho1Trans

  /**
   *  \brief Performs the transformation \f$ A' = O_1^T A O_1\f$ for
   *  a general matrix \f$A\f$ 
   *
   *  \f$ O_1 \f$ is the orthonormalization matrix stored in 
   *  AOIntegrals::ortho1 (see AOIntegrals::ORTHO_TYPE and
   *  AOIntegrals::computeOrtho for details).
   *
   */ 
  template <typename T> 
  void AOIntegrals::Ortho1TransT(T* A, T* TransA) {

    // Make sure that the incoming matricies are of the right size
    assert(memManager_.template getSize(A) == nSQ_);
    assert(memManager_.template getSize(TransA) == nSQ_);

    // Allocate scratch space
    T* SCR = memManager_.template malloc<T>(nSQ_);

    // If necessary, create a complex copy of ortho1
    T* O1 = reinterpret_cast<T*>(ortho1);
    if(std::is_same<T,dcomplex>::value) {
      O1 = memManager_.template malloc<T>(nSQ_);
      std::copy_n(ortho1,nSQ_,O1);
    }

    // Perform transformation
    Gemm('T', 'N', basisSet_.nBasis, basisSet_.nBasis, basisSet_.nBasis, T(1.),
      O1, basisSet_.nBasis, A, basisSet_.nBasis, T(0.), SCR, basisSet_.nBasis);
    Gemm('N', 'N', basisSet_.nBasis, basisSet_.nBasis, basisSet_.nBasis, T(1.),
      SCR, basisSet_.nBasis, O1, basisSet_.nBasis, T(0.), TransA,
      basisSet_.nBasis);

    // Free up scratch space
    memManager_.free(SCR);
    if(std::is_same<T,dcomplex>::value) memManager_.free(O1);
    
  }; // AOIntegrals::Ortho1TransT

}; // namespace ChronusQ

#endif
