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
    oper_t      fockScalar;
    oper_t      fockMz;
    oper_t      fockMy;
    oper_t      fockMx;
    oper_t_coll fock;

    // Orthonormal Fock
    oper_t      fockOrthoScalar;
    oper_t      fockOrthoMz;
    oper_t      fockOrthoMy;
    oper_t      fockOrthoMx;
    oper_t_coll fockOrtho;

    // Coulomb (J[P])
    oper_t JScalar;

    // Exchange (K[P])
    oper_t      KScalar;
    oper_t      KMz;
    oper_t      KMy;
    oper_t      KMx;
    oper_t_coll K;

    // Perturbation Tensor (G[P])
    oper_t      PTScalar;
    oper_t      PTMz;
    oper_t      PTMy;
    oper_t      PTMx;
    oper_t_coll PT;


    // Orthonormal density
    oper_t      onePDMOrthoScalar;
    oper_t      onePDMOrthoMz;
    oper_t      onePDMOrthoMy;
    oper_t      onePDMOrthoMx;
    oper_t_coll onePDMOrtho;




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
      PTScalar(nullptr), PTMz(nullptr), PTMy(nullptr), PTMx(nullptr),
      onePDMOrthoScalar(nullptr), onePDMOrthoMz(nullptr), 
        onePDMOrthoMy(nullptr), onePDMOrthoMx(nullptr) { 

      // Allocate SingleSlater Object
      alloc(); 

    }

    // See include/singleslater/impl.hpp for documentation 
    // on the following constructors

    // Same type
    SingleSlater(const SingleSlater &); // Copy
    SingleSlater(SingleSlater &&);      // Move

    // Different type
    template <typename U> SingleSlater(const SingleSlater<U> &); // Copy 
    template <typename U> SingleSlater(SingleSlater<U> &&);      // Move 

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
    void formDensity(){ };

    

  }; // class SingleSlater

}; // namespace ChronusQ

#endif
