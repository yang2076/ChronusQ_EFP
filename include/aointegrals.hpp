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


  enum CORE_HAMILTONIAN_TYPE {
    NON_RELATIVISTIC,
    EXACT_2C
  };


  /**
   *  The TwoBodyContraction struct. Stores information
   *  pertinant for a two body operator contraction with
   *  a one body (2 index) operator. z.B. The density matrix.
   */ 
  template <typename T, typename G>
  struct TwoBodyContraction {

    T*  X;  ///< 1-Body (2 index) operator to contraction
    G*  AX; ///< 1-Body (2 index) storage for the contraction

    bool HER; ///< Whether or not X is hermetian
    
    TWOBODY_CONTRACTION_TYPE contType;

  }; // struct TwoBodyContraction

  enum CONTRACTION_ALGORITHM {
    DIRECT,
    INCORE,
    DENFIT
  }; ///< 2-e Integral Contraction Algorithm

  enum ORTHO_TYPE {
    LOWDIN,
    CHOLESKY
  }; ///< Orthonormalization Scheme

  class AOIntegrals {
  public:

    typedef double* oper_t; ///< Storage of an operator
    typedef std::vector<oper_t> oper_t_coll; ///< A collection of operators

  private:


    size_t nTT_; ///< Reduced number of basis functions \f$ N_B(N_B+1)/2 \f$
    size_t nSQ_; ///< Squared basis functions \f$ N_B^2\f$
    size_t npTT_; ///< Reduced number of primitive functions \f$ N_P(N_P+1)/2 \f$
    size_t npSQ_; ///< Squared primitive functions \f$ N_P^2\f$

    CQMemManager &memManager_; ///< CQMemManager to allocate matricies
    Molecule     &molecule_;   ///< Molecule object for nuclear potential
    BasisSet     &basisSet_;   ///< BasisSet for the GTO basis defintion

    // General wrapper for 1-e integrals
    // See src/aointegrals/aointegrals_builders.cxx for documentation
    oper_t_coll OneEDriver(libint2::Operator, std::vector<libint2::Shell>&);

    // 1-e builder for in-house integral code
    template <size_t NOPER, bool SYMM, typename F>
    oper_t_coll OneEDriverLocal(const F&, std::vector<libint2::Shell>&);

    // local one body integrals

    // Overlap integrals
      
    // overlap integral of a shell pair  
    std::vector<std::vector<double>> computeOverlapS(libint2::ShellPair&,
                                         libint2::Shell&,libint2::Shell&);

    // horizontal recursion of contracted overlap integral 
    double hRRSab(libint2::ShellPair&, libint2::Shell&,libint2::Shell&,
                  int,int*,int,int*);

    // horizontal recursion of uncontracted overlap integral
    double hRRiPPSab(libint2::ShellPair::PrimPairData&,libint2::Shell&,libint2::Shell&,
                  int,int*,int,int*);

    // vertical recursion of uncontracted overlap integral
    double vRRSa0(libint2::ShellPair::PrimPairData&,libint2::Shell&,int,int*);

    // angular momentum integrals

    // angular momentum integrals of a shell pair
    std::vector<std::vector<double>> computeAngularL(libint2::ShellPair&,
                                        libint2::Shell&,libint2::Shell&);

    // vertical recursion of uncontracted angular momentum integral
    double Labmu(libint2::ShellPair::PrimPairData&,libint2::Shell&,libint2::Shell&,
                 double*,double*,int,int*,int,int*,int);

    // momentum integrals

    // electric dipole (velocity gauge) integrals of a shell pair
    std::vector<std::vector<double>> computeEDipoleE1_vel(libint2::ShellPair&,
                      libint2::Shell&, libint2::Shell&);

    // contracted momentum integral
    double Momentummu(libint2::ShellPair&,libint2::Shell&,libint2::Shell&,
                      int,int*,int,int*,int);

    // electric dipole integrals

    // electric dipole (length gauge) integrals of a shell pair
    std::vector<std::vector<double>> computeDipoleE1(libint2::ShellPair&,
                                         libint2::Shell&,libint2::Shell&);

    // contracted electric dipole integrals
    double DipoleE1(libint2::ShellPair&,libint2::Shell&,libint2::Shell&,
                    int,int*,int,int*,int);

    // electric quadrupole integrals

    // electric quadrupole integrals of a shell pair
    std::vector<std::vector<double>> computeEQuadrupoleE2_vel(libint2::ShellPair&,
                                                  libint2::Shell&,libint2::Shell&); 

    // contracted electric quadrupole integrals of a shell pair
    double QuadrupoleE2_vel( libint2::ShellPair&,libint2::Shell&,libint2::Shell&,
                             int,int*,int,int*,int,int );

    // magnetic dipole integrals

    // contracted magnetic dipole integral
    double MDipoleM1( libint2::ShellPair&,libint2::Shell&,libint2::Shell&,
                      int,int*,int,int*,int );
  
    // magnetic quadrupole integrals

    // contracted magnetic quadrupole integrals 
    double QuadrupoleM2_vel( libint2::ShellPair&,libint2::Shell&,libint2::Shell&,
                             int,int*,int,int*,int,int );

    // magnetic quadrupole integrals of a shell pair
    std::vector<std::vector<double>> computeMQuadrupoleM2_vel( libint2::ShellPair&,
                                                  libint2::Shell&,libint2::Shell&);

    // electric octupole integrals

    // electric octupole integrals of a shell pair
    std::vector<std::vector<double>> computeEOctupoleE3_vel( libint2::ShellPair&,
                                                 libint2::Shell&,libint2::Shell&);

    // contracted electric octupole integral
    double OctupoleE3_vel( libint2::ShellPair&,libint2::Shell&,libint2::Shell&,
                           int,int*,int,int*,int,int,int );


    // Taylor intrapolation of Boys function
    void computeFmTTaylor(double*,double,int,int);

    // nuclear potential integrals

    // contracted nuclear potential integrals of a shell pair
    std::vector<std::vector<double>> computePotentialV(
      const std::vector<libint2::Shell> &, libint2::ShellPair&, 
      libint2::Shell&,libint2::Shell&); 

    inline std::vector<std::vector<double>> computePotentialV(
      libint2::ShellPair& pair, libint2::Shell &s1, libint2::Shell &s2) {
    
      std::vector<libint2::Shell> dummy;
      return computePotentialV(dummy,pair,s1,s2);

    }

    // horizontal recursion of contracted nuclear potential integrals
    double hRRVab(const std::vector<libint2::Shell>&,libint2::ShellPair&,
                  libint2::Shell&,libint2::Shell&,int,int*,int,int*);

    // Bra vertical recursion of uncontracted nuclear potential integrals
    double vRRVa0(const std::vector<libint2::Shell>&,
                  libint2::ShellPair::PrimPairData&,libint2::Shell&,
                  double*,double*,int,int,int*,int);

    // horizontal recursion of uncontracted nuclear potential integrals
    double hRRiPPVab(const std::vector<libint2::Shell>&,
      libint2::ShellPair::PrimPairData&, libint2::Shell&,
      libint2::Shell&, int,int*,int,int*,double*,int,int);
    
    // Ket vertical recursion of uncontracted nuclear potential integrals
    double vRRV0b(const std::vector<libint2::Shell>&,
                  libint2::ShellPair::PrimPairData&,libint2::Shell&,
                  double*,double*,int,int,int*,int);

    // spin orbit integrals

    // spin orbit integrals of a shell pair
    std::vector<std::vector<double>> computeSL(
      const std::vector<libint2::Shell>&, libint2::ShellPair&,
      libint2::Shell&,libint2::Shell&);

    inline std::vector<std::vector<double>> computeSL(libint2::ShellPair &pair,
      libint2::Shell &s1, libint2::Shell &s2) {

      std::vector<libint2::Shell> dummy;
      return computeSL(dummy,pair,s1,s2);

    }

    // vertical recursion of uncontracted spin orbit integral
    double Slabmu(const std::vector<libint2::Shell>&, 
      libint2::ShellPair::PrimPairData&,libint2::Shell&,
      libint2::Shell&, double*,double*,int,int*,int,int*,int,int,int);

    // pV dot p integrals

    // pV dot p integrals of a shell pair
    std::vector<std::vector<double>> computepVdotp(
      const std::vector<libint2::Shell>&, libint2::ShellPair&,
      libint2::Shell&,libint2::Shell&);

    inline std::vector<std::vector<double>> computepVdotp(
      libint2::ShellPair &pair, libint2::Shell &s1, libint2::Shell &s2) {

      std::vector<libint2::Shell> dummy;
      return computepVdotp(dummy,pair,s1,s2);

    }

    // vertical recursion of uncontracted pV dot p integrals
    double pVpab(const std::vector<libint2::Shell>&,
      libint2::ShellPair::PrimPairData&,libint2::Shell&,
      libint2::Shell&, int,int*,int,int*,int,int); 

    // local one body integrals end

    public:

    // Control Variables
    CORE_HAMILTONIAN_TYPE coreType;
    CONTRACTION_ALGORITHM cAlg;      ///< Algorithm for 2-body contraction
    ORTHO_TYPE            orthoType; ///< Orthogonalization scheme

    double threshSchwartz; ///< Schwartz screening threshold



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

    oper_t_coll magDipole;     ///< Electric Dipole matrix     (length)
    oper_t_coll magQuadrupole; ///< Electric Quadrupole matrix (length)

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
      threshSchwartz(1e-12), cAlg(DIRECT), orthoType(LOWDIN), 
      memManager_(memManager), basisSet_(basis), molecule_(mol), 
      schwartz(nullptr), ortho1(nullptr), ortho2(nullptr), overlap(nullptr), 
      kinetic(nullptr), potential(nullptr), ERI(nullptr), coreType(NON_RELATIVISTIC) {

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


    // Print (see src/aointegrals/print.cxx for docs)
    friend std::ostream & operator<<(std::ostream &, const AOIntegrals& );


    // Memory

    // Deallocation (see src/aointegrals/aointegrals.cxx for docs)
    void dealloc();





    // Integral evaluation
    // (see src/aointegrals/aointegrals_builders.cxx for docs)

    void computeAOOneE(bool); // Evaluate the 1-e ints in the CGTO basis
    void computeERI();    // Evaluate and store the ERIs in the CGTO basis
    void computeOrtho();  // Evaluate orthonormalization transformations
    void computeSchwartz(); // Evaluate schwartz bounds over CGTOS

    // CH == Core Hamiltonian
    void computeCoreHam(CORE_HAMILTONIAN_TYPE); // Compute the CH
    void computeNRCH(double*); // Non-relativistic CH
    void computeX2CCH(std::vector<double*>&); // X2C CH (aointegrals_rel.cxx)
    void compute4CCH(std::vector<libint2::Shell>&, double *); // 4C CH

    // Allow for delayed evaluation of CH
    inline void computeCoreHam() { computeCoreHam(coreType); }

    // Integral contraction

    /**
     *  Contract the two body potential with one body (2 index) operators.
     *
     *  Smartly determines whether to do the contraction directly, incore
     *  or using density fitting depending on context
     *
     *  \param [in/ont] contList List of one body operators for contraction.
     */ 
    template <typename T, typename G>
    void twoBodyContract(std::vector<TwoBodyContraction<T,G>> &contList) {

      // Sanity check of dimensions
      assert( std::all_of(contList.begin(),contList.end(),
              [&](TwoBodyContraction<T,G> &C) {
                bool ret(true);
                ret = ret and (memManager_.template getSize<T>(C.X)  == nSQ_);
                ret = ret and (memManager_.template getSize<G>(C.AX) == nSQ_);
                return ret;
              }) );

      if( cAlg == INCORE ) twoBodyContractIncore(contList);
      else if( cAlg == DIRECT ) twoBodyContractDirect(contList);
    };
    

    // INCORE contraction routines
    // Perform the two body contraction incore (using the rank-4 ERI tensor)
    // see include/aointegrals/contract/incore.hpp for docs.
      
    template <typename T, typename G>
    void twoBodyContractIncore(std::vector<TwoBodyContraction<T,G>>&);

    template <typename T, typename G>
    void JContractIncore(TwoBodyContraction<T,G> &);

    template <typename T, typename G>
    void KContractIncore(TwoBodyContraction<T,G> &);



    // DIRECT contraction routines
    // Perform the two body contraction directly
    // see include/aointegrals/contract/direct.hpp for docs.
    template <typename T, typename G>
    void twoBodyContractDirect(std::vector<TwoBodyContraction<T,G>>&);

    template <typename T, typename G>
    void directScaffold(std::vector<TwoBodyContraction<T,G>>&);

    template <typename T, typename G>
    void JContractDirect(TwoBodyContraction<T,G> &);

    template <typename T, typename G>
    void KContractDirect(TwoBodyContraction<T,G> &);


    // Transformations to and from the orthonormal basis
    // see include/aointegrals/ortho.hpp for docs
      
    template <typename T> void Ortho1Trans(T* A, T* TransA); 
    template <typename T> void Ortho2Trans(T* A, T* TransA); 
    template <typename T> void Ortho1TransT(T* A, T* TransA);
    template <typename T> void Ortho2TransT(T* A, T* TransA);
    
  }; // class AOIntegrals

  extern std::array<std::array<double,25>,3201> FmTTable;

  void generateFmTTable();  


}; // namespace ChronusQ

#endif
