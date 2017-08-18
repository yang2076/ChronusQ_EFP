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
#ifndef __INCLUDED_SINGLESLATER_FOCK_HPP__
#define __INCLUDED_SINGLESLATER_FOCK_HPP__

#include <singleslater.hpp>
#include <cqlinalg/blasext.hpp>

namespace ChronusQ {

  /**
   *  \brief Forms the Fock matrix for a single slater determinant using
   *  the 1PDM.
   *
   *  \param [in] increment Whether or not the Fock matrix is being 
   *  incremented using a previous density
   *
   *  Populates / overwrites fock strorage
   */ 
  template <typename T>
  void SingleSlater<T>::formFock(bool increment) {

    size_t NB = aoints.basisSet().nBasis;
    size_t NB2 = NB*NB;

    // Form G[D]
    formGD(increment);

    // Zero out the Fock
    for(auto &F : fock) std::fill_n(F,NB2,0.);

    // Copy over the Core Hamiltonian
    SetMatRE('N',NB,NB,1.,aoints.coreH[SCALAR],NB,fock[SCALAR],NB);
    for(auto i = 1; i < aoints.coreH.size(); i++) 
      SetMatIM('N',NB,NB,1.,aoints.coreH[i],NB,fock[i],NB);

    // Add in the perturbation tensor
    for(auto i = 0ul; i < fock.size(); i++)
      MatAdd('N','N', NB, NB, T(1.), fock[i], NB, T(1.), GD[i], NB,
        fock[i], NB);

#if 0
    printFock(std::cout);
#endif

  }; // SingleSlater<T>::fockFock


  /**
   *  \brief Forms the Hartree-Fock perturbation tensor
   *
   *  Populates / overwrites GD storage (and JScalar and K storage)
   */ 
  template <typename T>
  void SingleSlater<T>::formGD(bool increment) {

    // Decide list of onePDMs to use
    oper_t_coll &contract1PDM  = increment ? deltaOnePDM : this->onePDM;


    size_t NB = aoints.basisSet().nBasis;
    size_t NB2 = NB*NB;

#if 0

    std::vector<TwoBodyContraction<T,double>> jContract =
      { {this->onePDM[SCALAR], JScalar, true, COULOMB} };
    

    std::vector<TwoBodyContraction<T,T>> kContract;

    // Determine how many (if any) exchange terms to calculate
    for(auto i = 0; i < K.size(); i++)
      kContract.push_back({this->onePDM[i], K[i], true, EXCHANGE});


    // Perform J contraction
    aoints.twoBodyContract(jContract);

    // Perform K contraction
    aoints.twoBodyContract(kContract);

#else




    // Possibly allocate a temporary for J matrix
    T* JContract;
    if(std::is_same<double,T>::value) 
      JContract = reinterpret_cast<T*>(JScalar);
    else {
      JContract = this->memManager.template malloc<T>(NB2);
    }

    // Zero out J
    if(not increment) memset(JContract,0,NB2*sizeof(T));

    std::vector<TwoBodyContraction<T,T>> contract =
      { {contract1PDM[SCALAR], JContract, true, COULOMB} };

    // Determine how many (if any) exchange terms to calculate
    for(auto i = 0; i < K.size(); i++) {
      contract.push_back({contract1PDM[i], K[i], true, EXCHANGE});

      // Zero out K[i]
      if(not increment) memset(K[i],0,NB2*sizeof(T));
    }

    aoints.twoBodyContract(contract);

    if(not std::is_same<double,T>::value) {
      if(not increment)
        GetMatRE('N',NB,NB,1.,JContract,NB,JScalar,NB);
      else {
        MatAdd('N','N',NB,NB,T(1.),JContract,NB,T(1.),
          JScalar,NB,JContract,NB);
        GetMatRE('N',NB,NB,1.,JContract,NB,JScalar,NB);
      }
      this->memManager.free(JContract);
    }

#endif

    // Form GD: G[D] = 2.0*J[D] - K[D]

    for(auto i = 0; i < K.size(); i++)
      MatAdd('N','N', NB, NB, T(0.), GD[i], NB, T(-1.), K[i], NB,
        GD[i], NB);

    MatAdd('N','N', NB, NB, T(1.), GD[SCALAR], NB, T(2.), JScalar, NB,
      GD[SCALAR], NB);
      
#if 0
    printJ(std::cout);
    printK(std::cout);
    printGD(std::cout);
#endif

  }; // SingleSlater<T>::formGD

}; // namespace ChronusQ

#endif
