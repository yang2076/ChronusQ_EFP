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
#include <util/matout.hpp>

namespace ChronusQ {

  /**
   *  \brief List of default atomic multiplicities for 
   *  atomic SCF.
   *
   *  XXX: Only supported through Kr
   */ 
  static std::map<size_t,size_t> defaultMult(
    {
      { 1  , 2 }, // H
      { 2  , 1 }, // He
      { 3  , 2 }, // Li
      { 4  , 1 }, // Be
      { 5  , 2 }, // B
      { 6  , 3 }, // C
      { 7  , 4 }, // N
      { 8  , 3 }, // O
      { 9  , 2 }, // F
      { 10 , 1 }, // Ne
      { 11 , 2 }, // Na
      { 12 , 1 }, // Mg
      { 13 , 2 }, // Al
      { 14 , 4 }, // Si
      { 15 , 6 }, // P
      { 16 , 3 }, // S
      { 17 , 2 }, // Cl
      { 18 , 1 }, // Ar
      { 19 , 2 }, // K
      { 20 , 1 }, // Ca
      { 21 , 2 }, // Sc
      { 22 , 3 }, // Ti
      { 23 , 4 }, // V
      { 24 , 7 }, // Cr
      { 25 , 2 }, // Mn
      { 26 , 5 }, // Fe
      { 27 , 4 }, // Co
      { 28 , 3 }, // Ni
      { 29 , 2 }, // Cu
      { 30 , 1 }, // Zn
      { 31 , 2 }, // Ga
      { 32 , 3 }, // Ge
      { 33 , 4 }, // As
      { 34 , 3 }, // Se
      { 35 , 2 }, // Br
      { 36 , 1 }, // Kr
    }
  );


  /**
   *  \brief Forms a set of guess orbitals for a single slater
   *  determininant SCF in various ways
   */ 
  template <typename T>
  void SingleSlater<T>::formGuess() {

    if( printLevel > 0 )
      std::cout << "  *** Forming Initial Guess Density for SCF Procedure ***"
                << std::endl << std::endl;

    if( aoints.molecule().nAtoms == 1  or scfControls.guess == CORE)
      CoreGuess();
    else if( scfControls.guess == SAD ) SADGuess();
    else if( scfControls.guess == RANDOM ) RandomGuess();

    // Common to all guess: form new set of orbitals from
    // initial guess at Fock.
    getNewOrbitals(false);
    
    // If RANDOM guess, scale the densites appropriately
    if( scfControls.guess == RANDOM ) {

      size_t NB    = aoints.basisSet().nBasis;

      double TS = 
        this->template computeOBProperty<double,SCALAR>(aoints.overlap);

      double TZ = 
        this->template computeOBProperty<double,MZ>(aoints.overlap);
      double TY = 
        this->template computeOBProperty<double,MY>(aoints.overlap);
      double TX = 
        this->template computeOBProperty<double,MX>(aoints.overlap);

      double magNorm = std::sqrt(TZ*TZ + TY*TY + TX*TX);

      Scale(NB*NB,T(this->nO)/TS,this->onePDM[SCALAR],1);

      for(auto k = 1; k < this->onePDM.size(); k++)
        if( magNorm > 1e-10 )
          Scale(NB*NB,T(this->nOA - this->nOB)/magNorm,this->onePDM[k],1);
        else
          std::fill_n(this->onePDM[k],NB*NB,0.);

        

    }

  }; // SingleSlater<T>::formGuess

  /**
   *  \brief Populates the initial Fock matrix with the core
   *  hamiltonian.
   *
   *  I.e. Neglecting electron-electron interaction. While this
   *  is the simplest choice for an initial guess, it is often
   *  a very poor guess.
   */
  template <typename T>
  void SingleSlater<T>::CoreGuess() {

    if( printLevel > 0 )
      std::cout << "    * Forming the Core Hamiltonian Guess (F = H)\n\n";

    size_t FSize = memManager.template getSize(fock[SCALAR]);
    size_t NB    = std::sqrt(FSize);

    // Zero out the Fock
    for(auto &F : fock) std::fill_n(F,FSize,0.);

    // Copy over the Core Hamiltonian
    SetMatRE('N',NB,NB,1.,aoints.coreH[SCALAR],NB,fock[SCALAR],NB);
    for(auto i = 1; i < aoints.coreH.size(); i++) 
      SetMatIM('N',NB,NB,1.,aoints.coreH[i],NB,fock[i],NB);

  }; // SingleSlater<T>::CoreGuess 

  /**
   *  
   */ 
  template <typename T>
  void SingleSlater<T>::SADGuess() {

    if( printLevel > 0 )
      std::cout << "    * Forming the Superposition of Atomic Densities Guess"
                << " (SAD)\n\n";

    size_t NB    = aoints.basisSet().nBasis;

    // Zero out the densities
    for(auto &X : this->onePDM) std::fill_n(X,NB*NB,0.);


    // Determine the unique atoms
    std::vector<Atom> uniqueElements(aoints.molecule().atoms);

    // Sort the Atoms by atomic number for std::unique
    std::sort(uniqueElements.begin(),uniqueElements.end(),
      [](Atom &a, Atom &b) {
        return a.atomicNumber < b.atomicNumber;
      });

    // Obtain unique elements on sorted list
    auto it = std::unique(uniqueElements.begin(),uniqueElements.end(),
                [](Atom &a, Atom &b){ 
                  return a.atomicNumber == b.atomicNumber;
                });

    // Remove excess elements
    uniqueElements.resize(std::distance(uniqueElements.begin(),it));

    if( printLevel > 0 )
      std::cout << "  *** Found " << uniqueElements.size() 
                << " unique Atoms in Molecule Specification ***\n";


    // Set all atomic centers to origin
    for(auto &u : uniqueElements) u.coord = {0.,0.,0.};

    // Create a map from atoms in molecule to unique atom index
    std::map<size_t,size_t> mapAtom2Uniq;
    for(auto iAtm = 0ul; iAtm < aoints.molecule().nAtoms; iAtm++) {
      size_t curAtomicNumber = aoints.molecule().atoms[iAtm].atomicNumber;

      auto el = std::find_if(uniqueElements.begin(),uniqueElements.end(),
                  [&](Atom &a){return a.atomicNumber == curAtomicNumber;});

      if( el == uniqueElements.end() ) 
        CErr("Whoops! Can't find atom in unique list! Contact a developer!");

      mapAtom2Uniq[iAtm] = std::distance(uniqueElements.begin(),el);
    }

    if( printLevel > 0 )
      std::cout << "  *** Running " << uniqueElements.size() 
                << " Atomic SCF calculations for SAD Guess ***\n\n";

    for(auto iUn = 0; iUn < uniqueElements.size(); iUn++) {

      // Get default multiplicity for the Atom
      size_t defaultMultip;
      try { 
        defaultMultip = defaultMult[uniqueElements[iUn].atomicNumber];
      } catch(...) {
        CErr("AtomZ = " + std::to_string(uniqueElements[iUn].atomicNumber) +
             " not supported for SAD Guess");
      }


      std::string multipName;
      switch( defaultMultip ) {

        case 1: multipName = "Singlet"; break;
        case 2: multipName = "Doublet"; break;
        case 3: multipName = "Triplet"; break;
        case 4: multipName = "Quadruplet"; break;
        case 5: multipName = "Quintuplet"; break;
        case 6: multipName = "Sextuplet"; break;
        case 7: multipName = "Septuplet"; break;
        case 8: multipName = "Octuplet"; break;
  
        default: multipName = "UNKNOWN"; break;

      }


      if( printLevel > 0 )
        std::cout << "    * Running AtomZ = " 
                  << uniqueElements[iUn].atomicNumber << " as a " 
                  << multipName << std::endl;

      Molecule atom(0,defaultMultip,{ uniqueElements[iUn] });
      BasisSet basis(aoints.basisSet().basisName, atom, 
                 aoints.basisSet().forceCart, false);
     
      AOIntegrals aointsAtom(this->memManager,atom,basis);
      
      aointsAtom.cAlg           = INCORE;
      aointsAtom.computeERI();
      aointsAtom.computeCoreHam();

      std::shared_ptr<SingleSlater<T>> ss;
      
      try {
        // Trigger error if not an HF object
        dynamic_cast<HartreeFock<T>&>(*this);
  
        ss = std::dynamic_pointer_cast<SingleSlater<T>> (
               std::make_shared<HartreeFock<T>>(
                 aointsAtom,1, ( defaultMultip == 1 )
               )
             );

        ss->printLevel = 0;
        ss->scfControls.doIncFock = false;        
        ss->scfControls.dampError = 1e-4;
        ss->scfControls.nKeep     = 8;

        ss->formGuess();
        ss->SCF();

      } catch( const std::bad_cast &e ) {

      }     

      size_t NBbasis = basis.nBasis;

      // Place the atomic densities into the guess density

      for(auto iAtm = 0; iAtm < mapAtom2Uniq.size(); iAtm++)
        if( mapAtom2Uniq[iAtm] == iUn ) {
          SetMat('N',NBbasis,NBbasis,T(1.),ss->onePDM[SCALAR],NBbasis,
            this->onePDM[SCALAR] + aoints.basisSet().mapCen2BfSt[iAtm]*(1+NB),
            NB);
        }

    }

    // Spin-Average the SAD density
    if( this->onePDM.size() > 1 )
      SetMat('N',NB,NB,T(this->nOA - this->nOB) / this->nO, 
        this->onePDM[SCALAR],NB,this->onePDM[MZ],NB);

    if( printLevel > 0 )
      std::cout << std::endl
                << "  *** Forming Initial Fock Matrix from SAD Density ***\n\n";

    formFock();

  }; // SingleSlater<T>::SADGuess



  template <typename T>
  void SingleSlater<T>::RandomGuess() {

    size_t NB    = aoints.basisSet().nBasis;

    // Set up random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
//  std::uniform_real_distribution<> dis(-5,5);
    std::normal_distribution<> dis(0,5);

    for(auto &F : this->fock) std::fill_n(F,NB*NB,0.);


    // Copy over the Core Hamiltonian
    SetMatRE('N',NB,NB,1.,aoints.coreH[SCALAR],NB,fock[SCALAR],NB);
    for(auto i = 1; i < aoints.coreH.size(); i++) 
      SetMatIM('N',NB,NB,1.,aoints.coreH[i],NB,fock[i],NB);
   

    // Randomize the Fock matricies
    for(auto &F : this->fock) {
      for(auto k = 0; k < NB*NB; k++) F[k] = dis(gen);
      HerMat('L',NB,F,NB);
    }

  }

}; // namespace ChronusQ

#endif
