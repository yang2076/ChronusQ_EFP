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
#ifndef __INCLUDED_AOINTEGRALS_HPP_
#define __INCLUDED_AOINTEGRALS_HPP_

#include <chronusq_sys.hpp>
#include <molecule.hpp>
#include <basisset/basisset_def.hpp>
#include <memmanager.hpp>
#include <libint2/engine.h>


namespace ChronusQ {

  class AOIntegrals {
  public:

    typedef double* oper_t; ///< Storage of an operator
    typedef std::vector<oper_t> oper_t_coll; ///< A collection of operators
  private:

    enum CONTRACTION_ALGORITHM {
      DIRECT,
      INCORE,
      DENFIT
    }; ///< 2-e Integral Contraction Algorithm

    enum ORTHO_TYPE {
      LOWDIN,
      CHOLESKY
    }; ///< Orthonormalization Scheme

    size_t nTT_; ///< Reduced number of basis functions \f$ N_B(N_B+1)/2 \f$
    size_t nSQ_; ///< Squared basis functions \f$ N_B^2\f$
    size_t npTT_; ///< Reduced number of primitive functions \f$ N_P(N_P+1)/2 \f$
    size_t npSQ_; ///< Squared primitive functions \f$ N_P^2\f$

    double threshSchwartz_; ///< Schwartz screening threshold

    CONTRACTION_ALGORITHM cAlg_;
    ORTHO_TYPE            orthoType_;

    CQMemManager &memManager_; ///< CQMemManager to allocate matricies
    Molecule     &molecule_;   ///< Molecule object for nuclear potential
    BasisSet     &basisSet_;   ///< BasisSet for the GTO basis defintion

    // General wrapper for 1-e integrals
    // See src/aointegrals/aointegrals_builders.cxx for documentation
    oper_t_coll OneEDriver(libint2::Operator, std::vector<libint2::Shell>&);

    public:

    enum ERI_CONTRACTION_TYPE {
      COULOMB, ///< (mn | kl) X(lk)
      EXCHANGE,///< (mn | kl) X(nk)
      PAIR     ///< (mn | kl) X(nl)
    }; ///< ERI Tensor Contraction Specification


    // Operator storage
      
    // Meta data relating to screening, orthonormalization, etc
      
    oper_t schwartz; ///< Schwartz bounds for the ERIs
    oper_t ortho1;
    oper_t ortho2;

    // 1-e integrals
    
    oper_t overlap;   ///< Overlap matrix 
    oper_t kinetic;   ///< Kinetic matrix 
    oper_t potential; ///< Nuclear potential matrix 

    oper_t_coll lenElecDipole;     ///< Electric Dipole matrix     (length)
    oper_t_coll lenElecQuadrupole; ///< Electric Quadrupole matrix (length)
    oper_t_coll lenElecOctupole;   ///< Electric Octuupole matrix  (length)

    oper_t_coll velElecDipole;     ///< Electric Dipole matrix     (velocity)
    oper_t_coll velElecQuadrupole; ///< Electric Quadrupole matrix (velocity)
    oper_t_coll velElecOctupole;   ///< Electric Octuupole matrix  (velocity)

    oper_t_coll coreH; ///< Core Hamiltonian (scalar and magnetization)
    
    
    // 2-e Storage
      
    oper_t ERI;     

    // Constructors
    
    // Disable default constructor
    AOIntegrals() = delete;

    AOIntegrals(CQMemManager &memManager, Molecule &mol, BasisSet &basis) :
      threshSchwartz_(1e-14), cAlg_(DIRECT), orthoType_(LOWDIN), 
      memManager_(memManager), basisSet_(basis), molecule_(mol), 
      schwartz(NULL), ortho1(NULL), ortho2(NULL), overlap(NULL), 
      kinetic(NULL), potential(NULL), ERI(NULL) {

      nTT_  = basis.nBasis * ( basis.nBasis + 1 ) / 2;
      nSQ_  = basis.nBasis * basis.nBasis;
      npTT_ = basis.nPrimitive * ( basis.nPrimitive + 1 ) / 2;
      npSQ_ = basis.nPrimitive * basis.nPrimitive;

    };

    // See src/aointegrals/aointegrals.cxx for documentation 
    // onf the following constructors

    AOIntegrals(const AOIntegrals &); // Copy constructor
    AOIntegrals(AOIntegrals &&);      // Move constructor

    /**
     *  Copy assignment.
     *
     *  Assign one AOIntegrals object to another by a copy
     */ 
    AOIntegrals& operator=(const AOIntegrals &) = default;
    
    /**
     *  Move assignment.
     *
     *  Assign one AOIntegrals object to another by a move
     */ 
    AOIntegrals& operator=(AOIntegrals &&) = default;


    /**
     *  Destructor.
     *
     *  Destructs an AOIntegrals object
     */ 
    ~AOIntegrals() { dealloc(); }
    

    // Member functions




    // Memory

    // Deallocation (see src/aointegrals/aointegrals.cxx for docs)
    void dealloc();





    // Integral evaluation / contraction
    // (see src/aointegrals/aointegrals_builders.cxx for docs)

    void computeAOOneE(); // Evaluate and store the 1-e ints in the CGTO basis
    void computeERI();    // Evaluate and store the ERIs in the CGTO basis
    void computeOrtho();  // Evaluate orthonormalization transformations

    
  }; // class AOIntegrals

}; // namespace ChronusQ

#endif
