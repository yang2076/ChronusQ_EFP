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
#ifndef __INCLUDED_SINGLESLATER_SCF_HPP__
#define __INCLUDED_SINGLESLATER_SCF_HPP__

#include <singleslater.hpp>
#include <util/matout.hpp>
#include <cqlinalg/blas1.hpp>

namespace ChronusQ {

  /**
   *  \brief Performs the self--consistant field procedure given a set of 
   *  orbitals.
   *
   *  \warning SCF procedure assumes that the 1PDM and orbital (mo1/2) storage
   *  has been populated in some way.
   */ 
  template <typename T>
  void SingleSlater<T>::SCF() {

    bool isConverged = false;
    computeEnergy();
    for( scfConv.nSCFIter = 0; scfConv.nSCFIter < scfControls.maxSCFIter; 
         scfConv.nSCFIter++) {

      // Save current state of the wave function (method specific)
      saveCurrentState();

      // Exit loop on convergence
      if(isConverged) break;

      // Form the Fock matrix D(k) -> F(k)
      formFock();

      // Get new orbtials and densities from current Fock: F(k) -> C/D(k + 1)
      getNewOrbitals();

      // Evaluate convergence
      isConverged = evalConver();

      // Print out iteration information
      printSCFProg(std::cout);

    }; // Iteration loop

  }; // SingleSlater<T>::SCF()



  /**
   *  \brief Print the current convergence information of the SCF
   *  procedure
   */ 
  template <typename T>
  void SingleSlater<T>::printSCFProg(std::ostream &out) {

    // SCF Iteration
    out << "  SCFIt: " <<std::setw(6) << std::left << scfConv.nSCFIter + 1;

    // Current Total Energy
    out << std::setw(18) << std::fixed << std::setprecision(10)
                         << std::left << this->totalEnergy;

    out << std::scientific << std::setprecision(7);
    // Current Change in Energy
    out << std::setw(14) << std::right << scfConv.deltaEnergy;
/*
    out << "   ";
    out << std::setw(13) << std::right << PARMS;
    if(!this->isClosedShell && this->nTCS_ == 1) {
      out << "   ";
      out << std::setw(13) << std::scientific << std::right 
                         << std::setprecision(7) << PBRMS;
    }
*/
  
    out << std::endl;
  }; // SingleSlater<T>::printSCFProg


  /**
   *  \brief Saves the current state of wave function
   *
   *  Saves a copy of the current AO 1PDM
   */ 
  template <typename T>
  void SingleSlater<T>::saveCurrentState() {
    size_t OSize = this->memManager.template getSize(fock[0]);

    // Copy over current AO density matrix
    for(auto i = 0; i < this->onePDM.size(); i++)
      std::copy_n(this->onePDM[i],OSize,curOnePDM[i]);

  }; // SingleSlater<T>::saveCurrentState()

  /**
   *  \brief Computes the change in the current wave function
   *
   *  Saves onePDM - curOnePDM in deltaOnePDM
   */ 
  template <typename T>
  void SingleSlater<T>::formDelta() {

    size_t NB = this->aoints.basisSet().nBasis;
    for(auto i = 0; i < this->onePDM.size(); i++)
      MatAdd('N','N',NB,NB,T(1.),this->onePDM[i],NB,T(-1.),
        curOnePDM[i],NB,deltaOnePDM[i],NB);

  }; // SingleSlater<T>:formDelta

  
  /**
   *  \brief Obtain a new set of orbitals given a Fock matrix.
   *
   *  Currently implements the fixed-point SCF procedure.
   */ 
  template <typename T>
  void SingleSlater<T>::getNewOrbitals() {

    // Transform AO fock into the orthonormal basis
    ao2orthoFock();

    // Diagonalize the orthonormal fock Matrix
    diagOrthoFock();

    // Form the orthonormal density (in the AO storage)
    formDensity();

    // Copy the AO storage to orthonormal storage and back transform
    // the density into the AO basis. This is because ortho2aoDen
    // requires the onePDMOrtho storage is populated.
    for(auto i = 0; i < this->onePDM.size(); i++)
      std::copy_n(this->onePDM[i],
        this->memManager.template getSize(onePDMOrtho[i]),
        onePDMOrtho[i]);

    // Transform the orthonormal density to the AO basis
    ortho2aoDen();

  }; // SingleSlater<T>::getNewOrbitals

  /**
   *  \brief Evaluate SCF convergence based on various criteria.
   *
   *  Checks the norm of [F,D], if converged -> SCF converged.
   *
   *  Checks change in energy and density between SCF iterations,
   *    if *both* converged -> SCF converged.
   */ 
  template <typename T>
  bool SingleSlater<T>::evalConver() {

    // Check energy convergence
      
    // Save copy of old Energy
    double oldEnergy = this->totalEnergy;

    // Compute new energy (with new Density)
    computeEnergy();
    scfConv.deltaEnergy = this->totalEnergy - oldEnergy;

    bool energyConv = std::abs(scfConv.deltaEnergy) < scfControls.eneConvTol;


    // Check density convergence

    formDelta(); // Get change in density
    double RMSScalar = 0.;

    bool denConv(true);

    // Check FP convergence
    bool FDConv(false);


    bool isConverged = FDConv or (energyConv and denConv);

    return isConverged;

  }; // SingleSlater<T>::evalConver


  /**
   *  \brief Diagonalize the orthonormal fock matrix
   *
   *  General purpose routine which diagonalizes the orthonormal
   *  fock matrix and stores a set of orthonormal MO coefficients
   *  (in WaveFunction::mo1 and possibly WaveFunction::mo2) and
   *  orbital energies. General for both 1 and 2 spin components
   */ 
  template <typename T>
  void SingleSlater<T>::diagOrthoFock() {

    size_t NB = this->aoints.basisSet().nBasis * this->nC;
    size_t NB2 = NB*NB;

    // Copy over the fockOrtho into MO storage
    if(this->nC == 1 and this->iCS) 
      std::transform(fockOrtho[0],fockOrtho[0] + NB2,this->mo1,
        [](T a){ return a / 2.; }
      );
    else if(this->nC == 1)
      for(auto j = 0; j < NB2; j++) {
        this->mo1[j] = 0.5 * (fockOrtho[0][j] + fockOrtho[1][j]); 
        this->mo2[j] = 0.5 * (fockOrtho[0][j] - fockOrtho[1][j]); 
      }
    else { 
      // TODO 2C 
    }

    // Diagonalize the Fock Matrix
    int INFO = HermetianEigen('V', 'L', NB, this->mo1, NB, this->eps1, 
      this->memManager );
    if( INFO != 0 ) CErr("HermetianEigen failed in Fock1",std::cout);

    if(this->nC == 1 and not this->iCS) {
      INFO = HermetianEigen('V', 'L', NB, this->mo2, NB, this->eps2, 
        this->memManager );
      if( INFO != 0 ) CErr("HermetianEigen failed in Fock2",std::cout);
    }

  }; // SingleSlater<T>::diagOrthoFock


  /**
   *  \brief Transforms all of the spin components of the AO fock
   *  matrix to the orthonormal basis.
   *
   *  Populates / overwrites fockOrtho storage
   */ 
  template <typename T>
  void SingleSlater<T>::ao2orthoFock() {

    for(auto i = 0; i < fock.size(); i++)
      this->aoints.Ortho1Trans(fock[i],fockOrtho[i]);

  }; // SingleSlater<T>::ao2orthoFock



  /**
   *  \brief Transforms all of the spin compoenents of the orthonormal
   *  1PDM to the AO basis.
   *
   *  Populates / overwrites onePDM storage
   */ 
  template <typename T>
  void SingleSlater<T>::ortho2aoDen() {

    for(auto i = 0; i < onePDMOrtho.size(); i++)
      this->aoints.Ortho1Trans(onePDMOrtho[i],this->onePDM[i]);

  }; // SingleSlater<T>::ao2orthoFock

}; // namespace ChronusQ

#endif
