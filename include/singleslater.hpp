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

// Debug print triggered by Wavefunction
  
#ifdef _WaveFunctionDebug
  #define _SingleSlaterDebug
#endif

namespace ChronusQ {

  /**
   *  \brief A struct to hold the information pertaining to
   *  the control of an SCF procedure.
   *
   *  Holds information like convergence critera, DIIS settings, 
   *  max iterations, etc.
   */ 
  struct SCFControls {

    // Convergence crieteria
    double denConvTol = 1e-8;  ///< Density convergence criteria
    double eneConvTol = 1e-10; ///< Energy convergence criteria

    // DIIS settings (TODO)

    // Misc control
    size_t maxSCFIter = 128; ///< Maximum SCF iterations.

  }; // SCFControls struct

  /**
   *  \brief A struct to hold the current status of an SCF procedure
   *
   *  Holds information like current density / energy changes, number of 
   *  iterations, etc.
   */ 
  struct SCFConvergence {

    double deltaEnergy;  ///< Convergence of Energy
    double RMSDenScalar; ///< RMS change in Scalar density
    double RMSDenMag;    ///< RMS change in magnetization (X,Y,Z) density

    size_t nSCFIter = 0; ///< Number of SCF Iterations

  }; // SCFConvergence struct

  template <typename T>
  class SingleSlater : public WaveFunction<T> {

  protected:
    // Useful typedefs
    typedef T*                   oper_t;
    typedef std::vector<oper_t>  oper_t_coll;

  private:
  public:
    
    // Operator storage

    // AO Fock Matrix
    oper_t      fockScalar; ///< Scalar AO Fock matrix
    oper_t      fockMz;     ///< Mz AO Fock matrix
    oper_t      fockMy;     ///< My AO Fock matrix
    oper_t      fockMx;     ///< Mx AO Fock matrix
    oper_t_coll fock;       ///< List of populated AO Fock matricies

    // Orthonormal Fock
    oper_t      fockOrthoScalar; ///< Scalar orthonormal Fock matrix
    oper_t      fockOrthoMz; ///< Mz orthonormal Fock matrix
    oper_t      fockOrthoMy; ///< My orthonormal Fock matrix
    oper_t      fockOrthoMx; ///< Mx orthonormal Fock matrix
    oper_t_coll fockOrtho;   ///< List of populated orthonormal Fock matricies

    // Coulomb (J[D])
    oper_t JScalar; ///< Scalar Coulomb Matrix

    // Exchange (K[D])
    oper_t      KScalar; ///< Scalar exact (HF) exchange matrix
    oper_t      KMz;     ///< Mz exact (HF) exchange matrix
    oper_t      KMy;     ///< My exact (HF) exchange matrix
    oper_t      KMx;     ///< Mx exact (HF) exchange matrix
    oper_t_coll K;       ///< List of populated exact (HF) exchange matricies

    // Exact Perturbation Tensor (G[D])
    oper_t      GDScalar;///< HF Scalar perturbation tensor
    oper_t      GDMz;    ///< HF Mz perturbation tensor
    oper_t      GDMy;    ///< HF My perturbation tensor
    oper_t      GDMx;    ///< HF Mx perturbation tensor
    oper_t_coll GD;      ///< List of populated HF perturbation tensors


    // Orthonormal density
    oper_t      onePDMOrthoScalar; ///< Scalar orthonormal 1PDM matrix
    oper_t      onePDMOrthoMz; ///< Mz orthonormal 1PDM matrix
    oper_t      onePDMOrthoMy; ///< My orthonormal 1PDM matrix
    oper_t      onePDMOrthoMx; ///< Mx orthonormal 1PDM matrix
    oper_t_coll onePDMOrtho;   ///< List of populated orthonormal 1PDM matricies


    // SCF Variables
    SCFControls    scfControls; ///< Controls for the SCF procedure
    SCFConvergence scfConv;     ///< Current status of SCF convergence


    // Constructors
      
    /**
     *  SingleSlater Constructor. Constructs a SingleSlater object
     *
     *  \param [in] aoi  AOIntegrals object (which handels the BasisSet, etc)
     *  \param [in] nC  Number of spin components (1 and 2 are supported)
     */ 
    SingleSlater(AOIntegrals &aoi, size_t nC) : 
      WaveFunction<T>(aoi,nC),
      fockScalar(nullptr), fockMz(nullptr), fockMy(nullptr), fockMx(nullptr),
      fockOrthoScalar(nullptr), fockOrthoMz(nullptr), fockOrthoMy(nullptr), 
        fockOrthoMx(nullptr),
      JScalar(nullptr),
      KScalar(nullptr), KMz(nullptr), KMy(nullptr), KMx(nullptr),
      GDScalar(nullptr), GDMz(nullptr), GDMy(nullptr), GDMx(nullptr),
      onePDMOrthoScalar(nullptr), onePDMOrthoMz(nullptr), 
        onePDMOrthoMy(nullptr), onePDMOrthoMx(nullptr) { 

      // Allocate SingleSlater Object
      alloc(); 

    }

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


    // Declarations from Quantum 
    // (see include/singleslater/quantum.hpp for docs)
    void formDensity();
    void computeEnergy();

    // Form a fock matrix (see include/singleslater/fock.hpp for docs)
    virtual void formFock(bool increment = false);
    void formGD();

    // Form initial guess orbitals
    // see include/singleslater/guess.hpp for docs)
    void formGuess();
    



    // SCF procedural functions (see include/singleslater/scf.hpp for docs)
      
    // Perform the SCF
    void SCF(); 

    // Transformation functions to and from the orthonormal basis
    void ao2orthoFock();
    void ortho2aoDen();

    // Evaluate convergence
    bool evalConver();

    // Obtain new orbitals
    void getNewOrbitals();

    // Misc procedural
    void diagOrthoFock();
    void printSCFProg(std::ostream &out = std::cout);

  }; // class SingleSlater

}; // namespace ChronusQ

#endif