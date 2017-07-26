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
#ifndef __INCLUDED_SINGLESLATER_QUANTUM_HPP__
#define __INCLUDED_SINGLESLATER_QUANTUM_HPP__

#include <singleslater.hpp>
#include <cqlinalg/blas3.hpp>
#include <quantum/properties.hpp>

namespace ChronusQ {

  /**
   *  \brief Forms the 1PDM using a set of orbitals 
   *
   *  specialization of Quantum::formDensity. Populates / overwrites
   *  onePDM storage
   */ 
  template <typename T>
  void SingleSlater<T>::formDensity() {

    size_t NB  = this->aoints.basisSet().nBasis * this->nC;
    size_t NB2 = NB*NB;

    if(this->nC == 1) { 

      // DS = DA = CA * CA**H
      Gemm('N', 'C', NB, NB, this->nOA, T(1.), this->mo1, NB, this->mo1,
        NB, T(0.), this->onePDM[SCALAR], NB);

      if(not this->iCS) {

        // DZ = DB = CB * CB**H
        Gemm('N', 'C', NB, NB, this->nOB, T(1.), this->mo2, NB, this->mo2,
          NB, T(0.), this->onePDM[MZ], NB);

        // DS = DA + DB
        // DZ = DA - DB
        for(auto j = 0; j < NB2; j++) {
          T tmp = this->onePDM[SCALAR][j];
          this->onePDM[SCALAR][j] = this->onePDM[SCALAR][j] + this->onePDM[MZ][j]; 
          this->onePDM[MZ][j]     = tmp - this->onePDM[MZ][j]; 
        }

      } else {

        // DS = 2 * DA
        std::transform(this->onePDM[SCALAR], this->onePDM[SCALAR] + NB2,
          this->onePDM[SCALAR],[](T a){ return 2.*a; }
        );

      }
    } else {
      // TODO: 2C
    }

  }; // SingleSlater<T>::formDensity


  /**
   *  \brief Computes the total energy of a single slater determinent
   *
   *  Given a 1PDM and a Fock matrix (specifically the core Hamiltonian and
   *  G[D]), compute the energy.
   *
   *  \warning Assumes that density and Fock matrix have the appropriate form
   *
   *  F(S) = F(A) + F(B)
   *  F(Z) = F(A) - F(B)
   *  ...
   *
   *  \f[
   *     E = \frac{1}{2} \mathrm{Tr}[\mathbf{P}(\mathbf{H} + \mathbf{F})] +
   *     V_{NN}
   *  \f]
   *
   *  Specialization of Quantum<T>::computeEnergy, populates / overwrites 
   *  OBEnergy and MBEnergy
   */ 
  template <typename T>
  void SingleSlater<T>::computeEnergy() {

    // Scalar core hamiltonian contribution to the energy
    this->OBEnergy = 
      this->template computeOBProperty<double,DENSITY_TYPE::SCALAR>(
        this->aoints.coreH[SCALAR]);

 

    // One body Spin Orbit
    dcomplex SOEnergy(0.,0.);
    if(this->aoints.coreH.size() > 1) {
      SOEnergy = 
        this->template computeOBProperty<dcomplex,DENSITY_TYPE::MZ>(
          this->aoints.coreH[MZ]);
      SOEnergy += 
        this->template computeOBProperty<dcomplex,DENSITY_TYPE::MY>(
          this->aoints.coreH[MY]);
      SOEnergy += 
        this->template computeOBProperty<dcomplex,DENSITY_TYPE::MX>(
          this->aoints.coreH[MX]);
    };

    SOEnergy *= dcomplex(0.,1.);
    this->OBEnergy -= std::real(SOEnergy);

    this->OBEnergy *= 0.5;


    // Compute many-body contribution to energy
    // *** These calls are safe as proper zeros are returned by
    // property engine ***
    this->MBEnergy = 
      this->template computeOBProperty<double,DENSITY_TYPE::SCALAR>(GD[SCALAR]);
    this->MBEnergy += 
      this->template computeOBProperty<double,DENSITY_TYPE::MZ>(GD[MZ]);
    this->MBEnergy += 
      this->template computeOBProperty<double,DENSITY_TYPE::MY>(GD[MY]);
    this->MBEnergy += 
      this->template computeOBProperty<double,DENSITY_TYPE::MX>(GD[MX]);

    this->MBEnergy *= 0.25;

    // Assemble total energy
    this->totalEnergy = 
      this->OBEnergy + this->MBEnergy + this->aoints.molecule().nucRepEnergy;

    // Sanity checks
    assert( not std::isnan(this->OBEnergy) );
    assert( not std::isnan(this->MBEnergy) );
    assert( not std::isinf(this->OBEnergy) );
    assert( not std::isinf(this->MBEnergy) );

  }; // SingleSlater<T>::computeEnergy

}; // namespace ChronusQ

#endif
