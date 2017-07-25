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
#ifndef __INCLUDED_AOINTEGRALS_HPP__
#define __INCLUDED_AOINTEGRALS_HPP__

#include <chronusq_sys.hpp>
#include <molecule.hpp>
#include <basisset/basisset_def.hpp>
#include <memmanager.hpp>
#include <libint2/engine.h>


namespace ChronusQ {

  /**
   *   \brief The AOIntegrals class. A class to handle the evaluation and
   *   storage of integrals in the atomic orbital (AO) gaussian--type
   *   orbital (GTO) basis. General functionality extends to both the 
   *   contracted (CGTO) and uncontracted (primitive, PGTO) bases.
   *
   *   Also handles the contraction of 2-body (3,4 index) integrals
   *   with 1-body (2 index) operators.
   *
   *   \warning Currently only functional for real GTOs.
   */

  enum TWOBODY_CONTRACTION_TYPE {
    COULOMB, ///< (mn | kl) X(lk)
    EXCHANGE,///< (mn | kl) X(nk)
    PAIR     ///< (mn | kl) X(nl)
  }; ///< 2-Body Tensor Contraction Specification


  /**
   *  The TwoBodyContraction struct. Stores information
   *  pertinant for a two body operator contraction with
   *  a one body (2 index) operator. z.B. The density matrix.
   */ 
  template <typename T>
  struct TwoBodyContraction {

    T*  X;  ///< 1-Body (2 index) operator to contraction
    T*  AX; ///< 1-Body (2 index) storage for the contraction

    bool HER; ///< Whether or not X is hermetian
    
    TWOBODY_CONTRACTION_TYPE contType;

  }; // struct TwoBodyContraction

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

    CONTRACTION_ALGORITHM cAlg_;      ///< Algorithm for 2-body contraction
    ORTHO_TYPE            orthoType_; ///< Orthogonalization scheme

    CQMemManager &memManager_; ///< CQMemManager to allocate matricies
    Molecule     &molecule_;   ///< Molecule object for nuclear potential
    BasisSet     &basisSet_;   ///< BasisSet for the GTO basis defintion

    // General wrapper for 1-e integrals
    // See src/aointegrals/aointegrals_builders.cxx for documentation
    oper_t_coll OneEDriver(libint2::Operator, std::vector<libint2::Shell>&);

    public:



    // Operator storage
      
    // Meta data relating to screening, orthonormalization, etc
      
    oper_t schwartz; ///< Schwartz bounds for the ERIs
    oper_t ortho1;   ///< Orthogonalization matrix which S -> I
    oper_t ortho2;   ///< Inverse of ortho1

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
      
    oper_t ERI;    ///< Electron-Electron repulsion integrals (4 index) 

    // Constructors
    
    // Disable default constructor
    AOIntegrals() = delete;

    /**
     *  AOIntegrals Constructor. Constructs an AOIntegrals object.
     *
     *  \param [in] memManager Memory manager for matrix allocation
     *  \param [in] mol        Molecule object for molecular specification
     *  \param [in] basis      The GTO basis for integral evaluation
     */ 
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


    // Getters
    CQMemManager& memManager() { return memManager_; }
    BasisSet&     basisSet()   { return basisSet_;   }
    Molecule&     molecule()   { return molecule_;   }


    // Memory

    // Deallocation (see src/aointegrals/aointegrals.cxx for docs)
    void dealloc();





    // Integral evaluation
    // (see src/aointegrals/aointegrals_builders.cxx for docs)

    void computeAOOneE(); // Evaluate and store the 1-e ints in the CGTO basis
    void computeERI();    // Evaluate and store the ERIs in the CGTO basis
    void computeOrtho();  // Evaluate orthonormalization transformations



    // Integral contraction

    /**
     *  Contract the two body potential with one body (2 index) operators.
     *
     *  Smartly determines whether to do the contraction directly, incore
     *  or using density fitting depending on context
     *
     *  \param [in/ont] contList List of one body operators for contraction.
     */ 
    template <typename T>
    void twoBodyContract(std::vector<TwoBodyContraction<T>> &contList) {
      twoBodyContractIncore(contList);
    };
    
    // Perform the two body contraction incore (using the rank-4 ERI tensor)
    // see include/aointegrals/contract.hpp for docs.
    template <typename T>
    void twoBodyContractIncore(std::vector<TwoBodyContraction<T>>&);



    // Transformations to and from the orthonormal basis
    // see include/aointegrals/ortho.hpp for docs
      
    template <typename T> void Ortho1Trans(T* A, T* TransA); 
    template <typename T> void Ortho2Trans(T* A, T* TransA); 
    template <typename T> void Ortho1TransT(T* A, T* TransA);
    template <typename T> void Ortho2TransT(T* A, T* TransA);
    
  }; // class AOIntegrals

}; // namespace ChronusQ

#endif
