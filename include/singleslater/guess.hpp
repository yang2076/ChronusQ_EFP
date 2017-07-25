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
#ifndef __INCLUDED_SINGLESLATER_GUESS_HPP__
#define __INCLUDED_SINGLESLATER_GUESS_HPP__

#include <singleslater.hpp>
#include <cqlinalg.hpp>

namespace ChronusQ {

  /**
   *  \brief Forms a set of guess orbitals for a single slater
   *  determininant in various ways
   *
   *  XXX: Currently only supports CORE guess
   */ 
  template <typename T>
  void SingleSlater<T>::formGuess() {

    size_t FSize = this->memManager.template getSize(fock[0]);

    // CORE guess (F = H)
      
    // Zero out the Fock
    for(auto &F : fock) std::fill_n(F,FSize,0.);

    // Copy over the Core Hamiltonian
    std::copy_n(this->aoints.coreH[0], FSize, fock[0]);
    // FIXME: This must be multiplied by "i" for 2C
    for(auto i = 1; i < this->aoints.coreH.size(); i++)
      std::copy_n(this->aoints.coreH[i], FSize, fock[i]);


    // Common to all guess
    getNewOrbitals();

  }; // SingleSlater<T>::formGuess

}; // namespace ChronusQ

#endif
