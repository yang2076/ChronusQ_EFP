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

    size_t NB = this->aoints.basisSet().nBasis * this->nC;
    size_t NB2 = NB*NB;

    // Form G[D]
    formGD();

    // Zero out the Fock
    for(auto &F : fock) std::fill_n(F,NB2,0.);

    // Copy over the Core Hamiltonian
    std::copy_n(this->aoints.coreH[SCALAR], NB2, fock[SCALAR]);
    // FIXME: This must be multiplied by "i" for 2C
    for(auto i = 1; i < this->aoints.coreH.size(); i++)
      std::copy_n(this->aoints.coreH[i], NB2, fock[i]);

    // Add in the perturbation tensor
    for(auto i = 0ul; i < fock.size(); i++)
      MatAdd('N','N', NB, NB, T(1.), fock[i], NB, T(1.), GD[i], NB,
        fock[i], NB);

  }; // SingleSlater<T>::fockFock


  /**
   *  XXX: formGD should be SingleSlater pure virtual and specialized
   *  in derived classes
   *
   *  \brief Forms the Hartree-Fock perturbation tensor
   *
   *  Populates / overwrites GD storage (and JScalar and K storage)
   */ 
  template <typename T>
  void SingleSlater<T>::formGD() {

    std::vector<TwoBodyContraction<T>> cont;
    
    // Always do Scalar Coulomb
    cont.push_back({this->onePDM[SCALAR], JScalar, true, COULOMB}); 

    // Determine how many (if any) exchange terms to calculate
    for(auto i = 0; i < K.size(); i++)
      cont.push_back({this->onePDM[i], K[i], true, EXCHANGE});

    // Perform contraction
    this->aoints.twoBodyContract(cont);

    // Form GD: G[D] = J[D] - 0.5*K[D]
    size_t NB = this->aoints.basisSet().nBasis * this->nC;
    size_t NB2 = NB*NB;

    for(auto i = 0; i < GD.size(); i++)
      MatAdd('N','N', NB, NB, T(0.), GD[i], NB, T(-1.), K[i], NB,
        GD[i], NB);

    MatAdd('N','N', NB, NB, T(1.), GD[SCALAR], NB, T(2.), JScalar, NB,
      GD[SCALAR], NB);
      

  }; // SingleSlater<T>::formGD

}; // namespace ChronusQ

#endif
