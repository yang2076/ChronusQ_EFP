/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2017 Li Research Group (University of Washington)
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
#ifndef __INCLUDED_WAVEFUNCTION_IMPL_HPP__
#define __INCLUDED_WAVEFUNCTION_IMPL_HPP__

#include <wavefunction.hpp>
#include <util/preprocessor.hpp>

// Template for a collective operation on the members of a 
// WaveFunction object
  
#define WaveFunction_COLLECTIVE_OP(OP_MEMBER,OP_OP) \
  /* Handle densities */\
  OP_OP(T,this,other,this->memManager,mo1); \
  OP_OP(T,this,other,this->memManager,mo2); \
  OP_OP(double,this,other,this->memManager,eps1); \
  OP_OP(double,this,other,this->memManager,eps2); \
  \
  /* Handle Member data */\
  OP_MEMBER(this,other,nO); OP_MEMBER(this,other,nV); \
  OP_MEMBER(this,other,nOA); OP_MEMBER(this,other,nOB);\
  OP_MEMBER(this,other,nVA); OP_MEMBER(this,other,nVB);\



namespace ChronusQ {

  /**
   *  Constructs a WaveFunction object from another of a another (possibly the 
   *  same) type by copy.
   *
   *  \param [in] other WaveFunction object to copy
   *  \param [in] dummy Dummy argument to fix calling signature for delegation 
   *    to copy constructor
   */ 
  template <typename T>
  template <typename U>
  WaveFunction<T>::WaveFunction(const WaveFunction<U> &other, int dummy) : 
    Quantum<T>(dynamic_cast<const Quantum<U>&>(other)),
    aoints(other.aoints) {

#ifdef _WaveFunctionDebug
    std::cout << "WaveFunction<T>::WaveFunction(const WaveFunction<U>&) "
              << "(this = " << this << ", other = " << &other << ")" 
              << std::endl;
#endif

    WaveFunction_COLLECTIVE_OP(COPY_OTHER_MEMBER,COPY_OTHER_MEMBER_OP);

  }; // WaveFunction<T>::WaveFunction(const WaveFunction<U> &)


  /**
   *  Constructs a WaveFunction object from another of a another (possibly the 
   *  same) type by move.
   *
   *  \warning Deallocates the passed WaveFunction object
   *
   *  \param [in] other WaveFunction object to move
   *  \param [in] dummy Dummy argument to fix calling signature for delegation 
   *    to move constructor
   */ 
  template <typename T>
  template <typename U>
  WaveFunction<T>::WaveFunction(WaveFunction<U> &&other, int dummy) : 
    Quantum<T>(dynamic_cast<Quantum<U>&&>(std::move(other))),
    aoints(other.aoints) {

#ifdef _WaveFunctionDebug
    std::cout << "WaveFunction<T>::WaveFunction(WaveFunction<U>&&) "
              << "(this = " << this << ", other = " << &other << ")" 
              << std::endl;
#endif

    WaveFunction_COLLECTIVE_OP(MOVE_OTHER_MEMBER,MOVE_OTHER_MEMBER_OP);

  }; // WaveFunction<T>::WaveFunction(WaveFunction<U> &&)


  // Delagate the copy constructor to the conversion constructors
  template <typename T>
  WaveFunction<T>::WaveFunction(const WaveFunction<T> &other) : 
    WaveFunction(other,0){ };
  template <typename T>
  WaveFunction<T>::WaveFunction(WaveFunction<T> &&other) : 
    WaveFunction(std::move(other),0){ };





  /**
   *  Allocates the internal memory a WaveFunction object
   */ 
  template <typename T>
  void WaveFunction<T>::alloc() {

#ifdef _WaveFunctionDebug
    std::cout << "WaveFunction::alloc (this = " << this << ")" << std::endl;
#endif

    size_t NB = this->nC * aoints.basisSet().nBasis;

    mo1  = this->memManager.template malloc<T>(NB*NB);
    eps1 = this->memManager.template malloc<double>(NB);

    if( this->nC == 1 and (not this->iCS) ) {
      mo2  = this->memManager.template malloc<T>(NB*NB);
      eps2 = this->memManager.template malloc<double>(NB);
    }

  }; // WaveFunction<T>::alloc



  /**
   *  Deallocates the internal memory a WaveFunction object
   */ 
  template <typename T>
  void WaveFunction<T>::dealloc() {

#ifdef _WaveFunctionDebug
    std::cout << "WaveFunction::dealloc (this = " << this << ")" << std::endl;
#endif

    WaveFunction_COLLECTIVE_OP(DUMMY3,DEALLOC_OP_5);

  }; // WaveFunction<T>::dealloc

}; // namespace ChronusQ
#endif
