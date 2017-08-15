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
#ifndef __INCLUDED_QUANTUM_PRINT_HPP__
#define __INCLUDED_QUANTUM_PRINT_HPP__

#include <quantum.hpp>
#include <util/matout.hpp>


namespace ChronusQ {

  template <typename T>
  void Quantum<T>::print1PDM(std::ostream &out) {

    size_t NB = std::sqrt(memManager.template getSize<T>(onePDM[0]));

    prettyPrintSmart(out,"1PDM (AO) Scalar",onePDM[SCALAR],NB,NB,NB);

    if( onePDM.size() > 1 )
      prettyPrintSmart(out,"1PDM (AO) MZ",onePDM[MZ],NB,NB,NB);

    if( onePDM.size() > 2 ) {
      prettyPrintSmart(out,"1PDM (AO) MY",onePDM[MY],NB,NB,NB);
      prettyPrintSmart(out,"1PDM (AO) MX",onePDM[MX],NB,NB,NB);
    }

  }; // Quantum<T>::print1PDM

};

#endif

