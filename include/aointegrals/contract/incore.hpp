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
#ifndef __INCLUDED_AOINTEGRALS_CONTRACT_INCORE_HPP__
#define __INCLUDED_AOINTEGRALS_CONTRACT_INCORE_HPP__


#include <aointegrals.hpp>
#include <util/threads.hpp>
#include <cqlinalg/blas3.hpp>

// Use stupid but bullet proof incore contraction for debug
//#define _BULLET_PROOF_INCORE

namespace ChronusQ {

  /**
   *  \brief Perform various tensor contractions of the full ERI
   *  tensor in core. Wraps other helper functions and provides
   *  loop structure
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
  template <typename T, typename G>
  void AOIntegrals::twoBodyContractIncore(
    std::vector<TwoBodyContraction<T,G>> &list) {

    auto topIncore = std::chrono::high_resolution_clock::now();

    // Loop over matricies to contract with
    for(auto &C : list) {

      // Coulomb-type (34,12) ERI contraction
      // AX(mn) = (mn | kl) X(kl)
      if( C.contType == COULOMB ) JContractIncore(C);

      // Exchange-type (23,12) ERI contraction
      // AX(mn) = (mk |ln) X(kl)
      else if( C.contType == EXCHANGE ) KContractIncore(C);

    } // loop over matricies

    auto botIncore = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> durIncore = botIncore - topIncore;
#ifdef _REPORT_INTEGRAL_TIMINGS
    std::cerr << "Incore Contraction took " << durIncore.count() << " s\n\n";
#endif
  }; // AOIntegrals::twoBodyContractIncore


  template<>
  void AOIntegrals::twoBodyContractIncore(
    std::vector<TwoBodyContraction<dcomplex,double>> &list) {

    assert(std::all_of(list.begin(),list.end(),
      [](TwoBodyContraction<dcomplex,double>& C){ 
        return C.contType == COULOMB; 
      }));

    for(auto &C : list) JContractIncore(C);
    
  }; // AOIntegrals::twoBodyContractIncore (complex, real)



  /**
   *  \brief Perform a Coulomb-type (34,12) ERI contraction with
   *  a one-body operator.
   */   
  template <typename T, typename G>
  void AOIntegrals::JContractIncore(TwoBodyContraction<T,G> &C) {

    // Temporaries Hermetian code
    double *X, *AX;
    if( C.HER ) {

      // Handle smart allocation of real matricies for Coulomb contraction
      // with Hermetian matricies
        
      if( std::is_same<double,T>::value ) 
        X = reinterpret_cast<double*>(C.X);
      else {
        // Copy over real part into allocated temporary
        X = memManager_.malloc<double>(nSQ_);
        std::transform(C.X,C.X + nSQ_,X,
          []( T a ) -> double { return std::real(a); }
        ); 
      }

      if( std::is_same<double,G>::value ) 
        AX = reinterpret_cast<double*>(C.AX);
      else {
        AX = memManager_.malloc<double>(nSQ_);  
        std::fill_n(AX,nSQ_,0.);
      }

    } else
      assert( (std::is_same<T,G>::value) );

    #ifdef _BULLET_PROOF_INCORE

    size_t NB3 = basisSet_.nBasis * nSQ_;

    // Hermetian code
    if(C.HER) 
    for(auto i = 0; i < basisSet_.nBasis; ++i)
    for(auto j = 0; j < basisSet_.nBasis; ++j)
    for(auto k = 0; k < basisSet_.nBasis; ++k)
    for(auto l = 0; l < basisSet_.nBasis; ++l) {
      AX[i + j*basisSet_.nBasis] +=
        ERI[i + j*basisSet_.nBasis + k*nSQ_ + l*NB3] *
        X[k + l*basisSet_.nBasis];
    }
    
    // Nonhermetian code
    else
    for(auto i = 0; i < basisSet_.nBasis; ++i)
    for(auto j = 0; j < basisSet_.nBasis; ++j)
    for(auto k = 0; k < basisSet_.nBasis; ++k)
    for(auto l = 0; l < basisSet_.nBasis; ++l) {
      reinterpret_cast<T*>(C.AX)[i + j*basisSet_.nBasis] +=
        ERI[i + j*basisSet_.nBasis + k*nSQ_ + l*NB3] *
        reinterpret_cast<T*>(C.X)[k + l*basisSet_.nBasis];
    }

    #else


    // Hermetian code
    if(C.HER) {

      Gemm('N','N',nSQ_,1,nSQ_,1.,ERI,nSQ_,X,nSQ_,0.,AX,nSQ_);

    // Non-hermetian code
    } else {

      Gemm('N','N',nSQ_,1,nSQ_,T(1.),ERI,nSQ_,reinterpret_cast<T*>(C.X),nSQ_,
        T(0.),reinterpret_cast<T*>(C.AX),nSQ_);

    }


    #endif

    // Cleanup temporaries
    if(C.HER) {
        
      if( not std::is_same<double,T>::value ) memManager_.free(X);
      if( not std::is_same<double,G>::value ) {
        // Copy over real result into persistant storage
        std::copy_n(AX,nSQ_,C.AX);
        memManager_.free(AX);
      }
    }

  }; // AOIntegrals::JContractIncore



  template <typename T, typename G>
  void AOIntegrals::KContractIncore(TwoBodyContraction<T,G> &C) {

    size_t NB3 = basisSet_.nBasis * nSQ_;

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

    size_t LAThreads = GetLAThreads();
    SetLAThreads(1);

    #pragma omp parallel for
    for(auto nu = 0; nu < basisSet_.nBasis; nu++) 
      Gemm('N','N',basisSet_.nBasis,1,nSQ_,T(1.),
        ERI  + nu * NB3,              basisSet_.nBasis,
        C.X,                          nSQ_,T(0.),
        C.AX + nu * basisSet_.nBasis, basisSet_.nBasis);

    SetLAThreads(LAThreads);

    #endif

  }; // AOIntegrals::KContractIncore

}; // namespace ChronusQ

#endif

