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
#ifndef __INCLUDED_SINGLESLATER_PRINT_HPP__
#define __INCLUDED_SINGLESLATER_PRINT_HPP__

#include <singleslater.hpp>
#include <util/matout.hpp>


namespace ChronusQ {

  template <typename T>
  void SingleSlater<T>::printFock(std::ostream &out) {

    size_t NB = aoints.basisSet().nBasis;

    prettyPrintSmart(out,"Fock (AO) Scalar",fock[SCALAR],NB,NB,NB);

    if( fock.size() > 1 )
      prettyPrintSmart(out,"Fock (AO) MZ",fock[MZ],NB,NB,NB);

    if( fock.size() > 2 ) {
      prettyPrintSmart(out,"Fock (AO) MY",fock[MY],NB,NB,NB);
      prettyPrintSmart(out,"Fock (AO) MX",fock[MX],NB,NB,NB);
    }


  }; // SingleSlater<T>::printFock

  template <typename T>
  void SingleSlater<T>::print1PDMOrtho(std::ostream &out) {

    size_t NB = aoints.basisSet().nBasis;

    prettyPrintSmart(out,"1PDM (Ortho) Scalar",onePDMOrtho[SCALAR],NB,NB,NB);

    if( onePDMOrtho.size() > 1 )
      prettyPrintSmart(out,"1PDM (Ortho) MZ",onePDMOrtho[MZ],NB,NB,NB);

    if( onePDMOrtho.size() > 2 ) {
      prettyPrintSmart(out,"1PDM (Ortho) MY",onePDMOrtho[MY],NB,NB,NB);
      prettyPrintSmart(out,"1PDM (Ortho) MX",onePDMOrtho[MX],NB,NB,NB);
    }


  }; // SingleSlater<T>::print1PDMOrtho

  template <typename T>
  void SingleSlater<T>::printGD(std::ostream &out) {

    size_t NB = aoints.basisSet().nBasis;

    prettyPrintSmart(out,"GD (AO) Scalar",GD[SCALAR],NB,NB,NB);

    if( GD.size() > 1 )
      prettyPrintSmart(out,"GD (AO) MZ",GD[MZ],NB,NB,NB);

    if( GD.size() > 2 ) {
      prettyPrintSmart(out,"GD (AO) MY",GD[MY],NB,NB,NB);
      prettyPrintSmart(out,"GD (AO) MX",GD[MX],NB,NB,NB);
    }


  }; // SingleSlater<T>::printGD


  template <typename T>
  void SingleSlater<T>::printJ(std::ostream &out) {

    size_t NB = aoints.basisSet().nBasis;

    prettyPrintSmart(out,"J (AO) Scalar",JScalar,NB,NB,NB);


  }; // SingleSlater<T>::printJ


  template <typename T>
  void SingleSlater<T>::printK(std::ostream &out) {

    size_t NB = aoints.basisSet().nBasis;

    prettyPrintSmart(out,"K (AO) Scalar",K[SCALAR],NB,NB,NB);

    if( K.size() > 1 )
      prettyPrintSmart(out,"K (AO) MZ",K[MZ],NB,NB,NB);

    if( K.size() > 2 ) {
      prettyPrintSmart(out,"K (AO) MY",K[MY],NB,NB,NB);
      prettyPrintSmart(out,"K (AO) MX",K[MX],NB,NB,NB);
    }


  }; // SingleSlater<T>::printK

}; // namespace ChronusQ

#endif
