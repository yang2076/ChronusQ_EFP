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

    printSCFHeader(std::cout);

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

    if(not isConverged)
      CErr(std::string("SCF Failed to converged within ") + 
        std::to_string(scfControls.maxSCFIter) + 
        std::string(" iterations"));
    else {
      std::cout << std::endl << "SCF Completed: E("
                << refShortName_ << ") = " << std::fixed
                << std::setprecision(10) << this->totalEnergy
                << " Eh after " << scfConv.nSCFIter
                << " SCF Iterations" << std::endl;
    } 
    std::cout << BannerEnd << std::endl;
    
  }; // SingleSlater<T>::SCF()

  
  template <typename T>
  void SingleSlater<T>::printSCFHeader(std::ostream &out) {
    out << BannerTop << std::endl;
    out << "Self Consistent Field (SCF) Settings:" << std::endl << std::endl;

    out << std::setw(38) << std::left << "  SCF Type:" << refLongName_ 
           << std::endl;

    out << std::setprecision(6) << std::scientific;
    out << std::setw(38)   << std::left << "  Density Convergence Tolerence:" 
           <<  scfControls.denConvTol << std::endl;

    out << std::setw(38)   << std::left << "  Energy Convergence Tolerence:" 
           <<  scfControls.eneConvTol << std::endl;

    out << std::setw(38) << std::left << "  Maximum Number of SCF Cycles:" 
           << scfControls.maxSCFIter << std::endl;

    out << std::endl << BannerMid << std::endl << std::endl;
    out << std::setw(16) << "SCF Iteration";
    out << std::setw(18) << "Energy (Eh)";
    out << std::setw(18) << "\u0394E (Eh)";
    out << std::setw(18) << " |\u0394P(S)|";
    if(not this->iCS or this->nC > 1)
      out << std::setw(18) << "  |\u0394P(M)|";
     
    out << std::endl;
    out << std::setw(16) << "-------------";
    out << std::setw(18) << "-----------";
    out << std::setw(18) << "-------";
    out << std::setw(18) << "-------";
    if(not this->iCS or this->nC > 1)
      out << std::setw(18) << "-------";
    out << std::endl;

  }; // SingleSlater<T>::printSCFHeader


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
    out << "   ";
    out << std::setw(13) << std::right << scfConv.RMSDenScalar;
    if(not this->iCS or this->nC > 1) {
      out << "   ";
      out << std::setw(13) << std::right << scfConv.RMSDenMag;
    }
  
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

    // Modify fock matrix if requested
    // TODO: need a better flag than this
    if (scfControls.doExtrap) modifyFock();

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

    size_t DSize = this->memManager. template getSize(fock[0]);
    scfConv.RMSDenScalar = TwoNorm<double>(DSize,deltaOnePDM[0],1);
    scfConv.RMSDenMag = 0.;
    for(auto i = 1; i < deltaOnePDM.size(); i++)
      scfConv.RMSDenMag += std::pow(TwoNorm<double>(DSize,deltaOnePDM[i],1),2.);
    scfConv.RMSDenMag = std::sqrt(scfConv.RMSDenMag);
    

    bool denConv = scfConv.RMSDenScalar < scfControls.denConvTol;

    // Check FP convergence
    bool FDConv(false);


    bool isConverged = FDConv or (energyConv and denConv);

    // Save matrices for next iteration if not converged
    if(not isConverged) saveSCFMatrices();

    // Toggle damping based on energy difference
    bool largeEDiff = std::abs(scfConv.deltaEnergy) > scfControls.dampError;
    if( scfControls.doDamp and not largeEDiff) {
      std::cout << "    *** Damping Disabled After "<<
        scfControls.dampError << " Converged Met ***" << std::endl;
      scfControls.dampParam = 0.;
    } else if( scfControls.doDamp and largeEDiff) {
      std::cout << "    *** Damping Enabled due to "<<
        scfControls.dampError << " Oscillation in Energy ***" << std::endl;
      scfControls.dampParam = scfControls.dampStartParam;
    }

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

 /**
   *  \brief Control routine for DIIS, damping and other
   *  ways to modify a Fock matrix during an SCF procedure
   *  
   */ 
  template <typename T>
  void SingleSlater<T>::modifyFock() {

    // DIIS extrapolation
//  if (scfControls.diisAlg != NONE) scfDIIS();

    // Static Damping
    if (scfControls.doDamp) fockDamping();

  }; // SingleSlater<T>::modifyFock

 /**
   *  \brief Static Fock Damping routine
   *
   *  F^k = (1-dp)*F^k + dp*F^{k-1}
   *  
   */ 
  template <typename T>
  void SingleSlater<T>::fockDamping() {

    std::cout << "  PJL inside fockDamping: " << std::endl;

    size_t NB = this->aoints.basisSet().nBasis;
    double dp = scfControls.dampParam;
   
    // Damp the current orthonormal fock matrix 
    for(auto i = 0; i < fockOrtho.size(); i++)
      MatAdd('N','N', NB, NB, T(1-dp), fockOrtho[i], NB, T(dp), 
        prevFock[i], NB, fockOrtho[i], NB);


  }; // SingleSlater<T>::fockDamping

 /**
   *  \brief Save SCF matrices for the next iteration
   *  
   */ 
  template <typename T>
  void SingleSlater<T>::saveSCFMatrices() {

    std::cout << "  PJL inside saveSCFMatrices: " << std::endl;

    // Copy the current orthonormal Fock to use during the next iteration 
    size_t NB = this->aoints.basisSet().nBasis;
    for(auto i = 0; i < this->fockOrtho.size(); i++)
      std::copy_n(this->fockOrtho[i],NB*NB,prevFock[i]);

  }; // SingleSlater<T>::saveSCFMatrices

}; // namespace ChronusQ

#endif
