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
#ifndef __INCLUDED_SINGLESLATER_HARTREEFOCK_HPP__
#define __INCLUDED_SINGLESLATER_HARTREEFOCK_HPP__

#include <chronusq_sys.hpp>
#include <singleslater.hpp>


namespace ChronusQ {


  template <typename T>
  class HartreeFock : public SingleSlater<T> {

    public:

    // Trivially inherit ctors from SingleSlater<T>

    template <typename... Args>
    HartreeFock(AOIntegrals &aoi, Args... args) : 
      SingleSlater<T>(aoi,args...), WaveFunctionBase(aoi,args...),
      QuantumBase(aoi.memManager(),args...) { 

      // Append HF tags to reference names
      if(this->nC == 1) {
        if(this->iCS) {
          this->refLongName_  += "Restricted Hartree-Fock";
          this->refShortName_ += "RHF";
        } else {
          this->refLongName_  += "Unrestricted Hartree-Fock";
          this->refShortName_ += "UHF";
        }
      } else {
        this->refLongName_  += "Generalized Hartree-Fock";
        this->refShortName_ += "GHF";
      }

    }; // HartreeFock constructor


    // Copy and Move ctors

    template <typename U>
    HartreeFock(const HartreeFock<U> &other, int dummy = 0) :
      SingleSlater<T>(dynamic_cast<const SingleSlater<U>&>(other)){ };
    template <typename U>
    HartreeFock(HartreeFock<U> &&other, int dummy = 0) :
      SingleSlater<T>(dynamic_cast<SingleSlater<U>&&>(std::move(other))){ };
    
    
    HartreeFock(const HartreeFock<T> &other, int dummy = 0) :
      SingleSlater<T>(dynamic_cast<const SingleSlater<T>&>(other)){ };
    HartreeFock(HartreeFock<T> &&other, int dummy = 0) :
      SingleSlater<T>(dynamic_cast<SingleSlater<T>&&>(std::move(other))){ };

  }; // class HartreeFock

}; // namespace ChronusQ

#endif
