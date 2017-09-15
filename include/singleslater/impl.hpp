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
#ifndef __INCLUDED_SINGLESLATER_IMPL_HPP__
#define __INCLUDED_SINGLESLATER_IMPL_HPP__

#include <singleslater.hpp>
#include <util/preprocessor.hpp>
#include <quantum/preprocessor.hpp>

// Template for a collective operation on the members of a 
// SingleSlater object
  

// FIXME: For copy and move, this only populates the lists, not the
// explicit pointers
#define SingleSlater_COLLECTIVE_OP(OP_OP,OP_VEC_OP) \
  /* Handle Operators */\
  OP_VEC_OP(T,this,other,memManager,fock); \
  OP_VEC_OP(T,this,other,memManager,fockOrtho); \
  OP_OP(double,this,other,memManager,JScalar); \
  OP_VEC_OP(T,this,other,memManager,K); \
  OP_VEC_OP(T,this,other,memManager,GD);\
  OP_VEC_OP(T,this,other,memManager,onePDMOrtho);\
  OP_VEC_OP(T,this,other,memManager,curOnePDM);\
  OP_VEC_OP(T,this,other,memManager,deltaOnePDM);

namespace ChronusQ {

  /**
   *  Constructs a SingleSlater object from another of a another (possibly the 
   *  same) type by copy.
   *
   *  \param [in] other SingleSlater object to copy
   *  \param [in] dummy Dummy argument to fix calling signature for delegation 
   *    to copy constructor
   */ 
  template <typename T>
  template <typename U>
  SingleSlater<T>::SingleSlater(const SingleSlater<U> &other,int dummy) : 
    QuantumBase(dynamic_cast<const QuantumBase&>(other)),
    WaveFunctionBase(dynamic_cast<const WaveFunctionBase&>(other)),
    SingleSlaterBase(dynamic_cast<const SingleSlaterBase&>(other)),
    WaveFunction<T>(dynamic_cast<const WaveFunction<U>&>(other)) {

#ifdef _SingleSlaterDebug
    std::cout << "SingleSlater<T>::SingleSlater(const SingleSlater<U>&) "
              << "(this = " << this << ", other = " << &other << ")" 
              << std::endl;
#endif

    SingleSlater_COLLECTIVE_OP(COPY_OTHER_MEMBER_OP, COPY_OTHER_MEMBER_VEC_OP);

  }; // SingleSlater<T>::SingleSlater(const SingleSlater<U> &)



  /**
   *  Constructs a SingleSlater object from another of a another (possibly the 
   *  same) by move.
   *
   *  \warning Deallocates the passed SingleSlater object
   *
   *  \param [in] other SingleSlater object to move
   *  \param [in] dummy Dummy argument to fix calling signature for delegation 
   *    to move constructor
   */ 
  template <typename T>
  template <typename U>
  SingleSlater<T>::SingleSlater(SingleSlater<U> &&other,int dummy) : 
    QuantumBase(dynamic_cast<QuantumBase&&>(std::move(other))),
    WaveFunctionBase(dynamic_cast<WaveFunctionBase&&>(std::move(other))),
    SingleSlaterBase(dynamic_cast<SingleSlaterBase&&>(std::move(other))),
    WaveFunction<T>(dynamic_cast<WaveFunction<U>&&>(std::move(other))) {

#ifdef _SingleSlaterDebug
    std::cout << "SingleSlater<T>::SingleSlater(SingleSlater<U>&&) "
              << "(this = " << this << ", other = " << &other << ")" 
              << std::endl;
#endif

    SingleSlater_COLLECTIVE_OP(MOVE_OTHER_MEMBER_OP, MOVE_OTHER_MEMBER_VEC_OP);

  }; // SingleSlater<T>::SingleSlater(SingleSlater<U> &&)


  // Delagate the copy constructor to the conversion constructors
  template <typename T>
  SingleSlater<T>::SingleSlater(const SingleSlater<T> &other) : 
    SingleSlater(other,0){ };
  template <typename T>
  SingleSlater<T>::SingleSlater(SingleSlater<T> &&other) : 
    SingleSlater(std::move(other),0){ };




  /**
   *  Allocates the internal memory a SingleSlater object
   */ 
  template <typename T>
  void SingleSlater<T>::alloc() {

#ifdef _SingleSlaterDebug
    std::cout << "SingleSlater::alloc (this = " << this << ")" << std::endl;
#endif

    size_t NB = aoints.basisSet().nBasis;

    SPIN_OPERATOR_ALLOC(NB,fock);
    SPIN_OPERATOR_ALLOC(NB,fockOrtho);
    SPIN_OPERATOR_ALLOC(NB,K);
    SPIN_OPERATOR_ALLOC(NB,GD);
    SPIN_OPERATOR_ALLOC(NB,onePDMOrtho);

    // J only has a scalar component
    JScalar = memManager.template malloc<double>(NB*NB);

    SPIN_OPERATOR_ALLOC(NB,curOnePDM);
    SPIN_OPERATOR_ALLOC(NB,deltaOnePDM);

  }; // SingleSlater<T>::alloc


  /**
   *  Deallocates the internal memory a SingleSlater object
   */ 
  template <typename T>
  void SingleSlater<T>::dealloc() {

#ifdef _SingleSlaterDebug
    std::cout << "SingleSlater::dealloc (this = " << this << ")" << std::endl;
#endif

    SingleSlater_COLLECTIVE_OP(DEALLOC_OP_5, DEALLOC_VEC_OP_5);

  }; // SingleSlater<T>::dealloc

}; // namespace ChronusQ


// Other implementation files
#include <singleslater/quantum.hpp> // Quantum declarations
#include <singleslater/fock.hpp>    // Fock matrix header
#include <singleslater/guess.hpp>   // Guess header
#include <singleslater/scf.hpp>     // SCF header
#include <singleslater/extrap.hpp>  // Extrapolate header
#include <singleslater/print.hpp>   // Print header
#include <singleslater/pop.hpp>     // Population analysis

#include <singleslater/kohnsham/impl.hpp> // KS headers

#endif
