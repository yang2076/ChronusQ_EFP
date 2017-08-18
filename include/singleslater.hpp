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
#ifndef __INCLUDED_SINGLESLATER_HPP__
#define __INCLUDED_SINGLESLATER_HPP__

#include <chronusq_sys.hpp>
#include <wavefunction.hpp>
#include <singleslater/base.hpp>

// Debug print triggered by Wavefunction
  
#ifdef _WaveFunctionDebug
  #define _SingleSlaterDebug
#endif

namespace ChronusQ {


  /**
   *  \brief The SingleSlater class. The typed abstract interface for all
   *  classes for which the wave function is described by a single slater
   *  determinant (HF, KS, PHF, etc).
   *
   *  Adds knowledge of storage type to SingleSlaterBase
   *
   *  Specializes the WaveFunction class of the same type
   */ 
  template <typename T>
  class SingleSlater : public SingleSlaterBase, public WaveFunction<T> {

  protected:

    // Useful typedefs
    typedef T*                        oper_t;
    typedef std::vector<oper_t>       oper_t_coll;
    typedef std::vector<oper_t_coll>  oper_t_coll2;

  private:
  public:
    
    // Operator storage

    // AO Fock Matrix
    oper_t_coll fock;       ///< List of populated AO Fock matricies

    // Orthonormal Fock
    oper_t_coll fockOrtho;   ///< List of populated orthonormal Fock matricies

    // Coulomb (J[D])
    double* JScalar; ///< Scalar Coulomb Matrix

    // Exchange (K[D])
    oper_t_coll K;       ///< List of populated exact (HF) exchange matricies

    // Exact Perturbation Tensor (G[D])
    oper_t_coll GD;      ///< List of populated HF perturbation tensors


    // Orthonormal density
    oper_t_coll onePDMOrtho;   ///< List of populated orthonormal 1PDM matricies


    // Current / change in state information (for use with SCF)
    oper_t_coll curOnePDM;    ///< List of the current 1PDMs
    oper_t_coll deltaOnePDM;  ///< List of the changes in the 1PDMs

    // Stores the previous Fock matrix to use for damping    
    oper_t_coll prevFock;     ///< AO Fock from the previous SCF iteration

    // Stores the previous matrices to use for DIIS 
    oper_t_coll2 diisFock;    ///< List of AO Fock matrices for DIIS extrap
    oper_t_coll2 diisOnePDM;  ///< List of AO Density matrices for DIIS extrap
    oper_t_coll2 diisError;   ///< List of orthonormal [F,D] for DIIS extrap

    // Constructors
      
    /**
     *  SingleSlater Constructor. Constructs a SingleSlater object
     *
     *  \param [in] aoi  AOIntegrals object (which handels the BasisSet, etc)
     *  \param [in] args Parameter pack for the remaining parameters of the
     *                   WaveFunction constructor. See include/wavefunction.hpp
     *                   for details. 
     */ 
    template <typename... Args>
    SingleSlater(AOIntegrals &aoi, Args... args) : 
      SingleSlaterBase(aoi,args...), WaveFunctionBase(aoi,args...),
      QuantumBase(aoi.memManager(),args...), WaveFunction<T>(aoi,args...), 
      JScalar(nullptr) {

      // Allocate SingleSlater Object
      alloc(); 

      // Determine Real/Complex part of method string
      if(std::is_same<T,double>::value) {
        refLongName_  = "Real ";
        refShortName_ = "R-";
      } else {
        refLongName_  = "Complex ";
        refShortName_ = "C-";
      }

    }; // SingleSlater constructor

    // See include/singleslater/impl.hpp for documentation 
    // on the following constructors

    // Different type
    template <typename U> SingleSlater(const SingleSlater<U> &, int dummy = 0);
    template <typename U> SingleSlater(SingleSlater<U> &&     , int dummy = 0);

    // Same type
    SingleSlater(const SingleSlater &);
    SingleSlater(SingleSlater &&);     

    /**
     *  Destructor.
     *
     *  Destructs a SingleSlater object
     */ 
    ~SingleSlater() { dealloc(); }



    // Public Member functions

    // Deallocation (see include/singleslater/impl.hpp for docs)
    void alloc();
    void dealloc();


    // Declarations from QuantumBase 
    // (see include/singleslater/quantum.hpp for docs)
    void formDensity();
    void computeEnergy();

    // Form a fock matrix (see include/singleslater/fock.hpp for docs)
    virtual void formFock(bool increment = false);
    void formGD(bool increment = false);

    // Form initial guess orbitals
    // see include/singleslater/guess.hpp for docs)
    void formGuess();
    



    // SCF procedural functions (see include/singleslater/scf.hpp for docs)

    // Transformation functions to and from the orthonormal basis
    void ao2orthoFock();
    void ortho2aoDen();

    // Evaluate convergence
    bool evalConver();

    // Obtain new orbitals
    void getNewOrbitals(bool frmFock = true);

    // Misc procedural
    void diagOrthoFock();
    void FDCommutator(oper_t_coll &);
    virtual void saveCurrentState();
    virtual void formDelta();
    void SCFInit();
    void SCFFin();

    // Print functions
    void printFock(std::ostream& )    ;
    void print1PDMOrtho(std::ostream&);
    void printGD(std::ostream&)       ;
    void printJ(std::ostream&)        ;
    void printK(std::ostream&)        ;

    // SCF extrapolation functions (see include/singleslater/extrap.hpp for docs)
    void allocExtrapStorage();
    void deallocExtrapStorage();
    void modifyFock();
    void fockDamping();
    void scfDIIS(size_t);

  }; // class SingleSlater

}; // namespace ChronusQ


// Include headers for specializations of SingleSlater
#include <singleslater/hartreefock.hpp> // HF specialization

#endif
