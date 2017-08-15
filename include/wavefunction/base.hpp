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
#ifndef __INCLUDED_WAVEFUNCTION_BASE_HPP__
#define __INCLUDED_WAVEFUNCTION_BASE_HPP__

#include <chronusq_sys.hpp>
#include <aointegrals.hpp>

namespace ChronusQ {

  /**
   *  \brief The WaveFunctionBase class. The abstraction of information
   *  relating to the WaveFunction class which are independent of storage
   *  type.
   *
   *  Specializes QuantumBase interface.
   *
   *  See WaveFunction for further docs.
   */ 
  class WaveFunctionBase : virtual public QuantumBase {

  public:
    // Member data

    AOIntegrals &aoints; ///< AOIntegrals for the evaluation of GTO integrals

    size_t nO;  ///< Total number of occupied orbitals
    size_t nV;  ///< Total number of virtual orbitals
    size_t nOA; ///< Number of occupied alpha orbitals (nC == 1)
    size_t nOB; ///< Number of occupied beta orbitals  (nC == 1)
    size_t nVA; ///< Number of virtual alpha orbitals  (nC == 1)
    size_t nVB; ///< Number of virtual beta orbitals   (nC == 1)

    
    // Disable default constructor
    WaveFunctionBase() = delete;

    // Default copy and move constructors
    WaveFunctionBase(const WaveFunctionBase &) = default;
    WaveFunctionBase(WaveFunctionBase &&)      = default;



    /**
     *  WaveFunctionBase Constructor. Constructs a WaveFunctionBase object
     *
     *  \param [in] aoi  AOIntegrals object (which handels the BasisSet, etc)
     *  \param [in] _nC  Number of spin components (1 and 2 are supported)
     *  \param [in] iCS  Whether or not to treat as closed shell
     */ 
    WaveFunctionBase(AOIntegrals &aoi, size_t _nC, bool iCS) : 
      QuantumBase(aoi.memManager(),_nC,iCS), aoints(aoi) {

      // Compute meta data

      nO = aoints.molecule().nTotalE;
      nV = 2*aoints.basisSet().nBasis - nO;

      if( iCS ) {
        nOA = nO / 2; nOB = nO / 2;
        nVA = nV / 2; nVB = nV / 2;
      } else {
        size_t nSingleE = aoints.molecule().multip - 1;
        nOB = (nO - nSingleE) / 2;
        nOA = nOB + nSingleE;
        nVA = aoints.basisSet().nBasis - nOA;
        nVB = aoints.basisSet().nBasis - nOB;
      }

    }; // WaveFunctionBase ctor


    // Print Functions
    virtual void printMO(std::ostream&)  = 0;
    virtual void printEPS(std::ostream&) = 0;

  }; // class WaveFunctionBase

}; // namespace ChronusQ

#endif
