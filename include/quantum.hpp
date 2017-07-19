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
#ifndef __INCLUDED_QUANTUM_HPP__
#define __INCLUDED_QUANTUM_HPP__

#include <chronusq_sys.hpp>
#include <memmanager.hpp>

// Debug print (triggers WaveFunction, etc)
#define _QuantumDebug

namespace ChronusQ {

  enum DENSITY_TYPE {
    SCALAR,MZ,MX,MY
  }; ///< Enumerate the types of densities for contraction

  template <typename T>
  class Quantum {
  protected:
    // Useful typedefs
    typedef T*                   oper_t;
    typedef std::vector<oper_t>  oper_t_coll;

    typedef std::array<double,3>    cartvec_t;
    typedef std::array<cartvec_t,3> cartmat_t;
    typedef std::array<cartmat_t,3> cartrk3_t;
    

  private:

    // Functions for the automatic evaluation of properties
    // see include/quantum/properties.hpp for documentation

    template<typename Scalar, typename Left, typename Right>
    static Scalar OperatorTrace(const Left& , const Right& ); 

    template<typename Scalar, DENSITY_TYPE DenTyp, typename Op>
    Scalar OperatorSpinCombine(const Op&);

  public:

    CQMemManager& memManager; ///< Memory manager for matrix allocation

    int   nC;   ///< Number of spin components
    bool  iCS;  ///< is closed shell?

    // 1PDM storage

    oper_t onePDMScalar; ///< Scalar part of the 1PDM
    oper_t onePDMMz;     ///< Z-Magnetization part of the 1PDM
    oper_t onePDMMy;     ///< Y-Magnetization part of the 1PDM
    oper_t onePDMMx;     ///< X-Magnetization part of the 1PDM
    
    oper_t_coll onePDM;  ///< 1PDM array (Scalar + Magnetization)

    // Property storage

    // Length gauge electric multipoles
    cartvec_t elecDipole;     ///< Electric Dipole in the length gauge
    cartmat_t elecQuadrupole; ///< Electric Quadrupole in the length gauge
    cartrk3_t elecOctupole;   ///< Electric Octupole in the length gauge
    
    // Spin expectation values
    cartvec_t SExpect; ///< Expectation values of Sx, Sy and Sz
    double    SSq;     ///< Expectation value of S^2


    // Constructors
      
    // Diable default constructor
    Quantum() = delete;

    /**
     *  Quantum Constructor. Constructs a Quantum object.
     *
     *  \param [in] mem   CQMemManager to handle to allocation of densities
     *  \param [in] _nC   Number of spin components (1 and 2 are supported)
     *  \param [in] _iCS  Whether or not system is closed shell
     *                    (only used when _nC == 1)
     *  \param [in] N     Dimension of the density matricies to be allocated
     */ 
    Quantum(CQMemManager &mem, size_t _nC = 1, bool _iCS = true, 
      size_t N = 0): 
      memManager(mem), nC(_nC), iCS(_iCS), onePDMScalar(nullptr), 
      onePDMMz(nullptr), onePDMMx(nullptr), onePDMMy(nullptr), 
      elecDipole({0.,0.,0.}),
      elecQuadrupole{{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
      elecOctupole{
        {
          {{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
          {{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
          {{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}}
        }}, SExpect({0.,0.,0.}), SSq(0.) {

        // Allocate densities
        if( N != 0 ) alloc(N);

    }; // Quantum::Quantum
    
       
    // See include/quantum/impl.hpp for documentation 
    // on the following constructors

    // Same type
    Quantum(const Quantum &); // Copy constructor
    Quantum(Quantum &&);      // Move constructor

    // Different type
    template <typename U> Quantum(const Quantum<U> &); // Copy constructor
    template <typename U> Quantum(Quantum<U> &&);      // Move constructor

    /**
     *  Deconstructor
     */ 
    ~Quantum(){ dealloc(); }

    // Member functions

    // Deallocation (see include/quantum/impl.hpp for docs)
    void alloc(size_t);
    void dealloc();




    // Procedural
      
    /**
     *  Function to form the density. Pure virtual
     */ 
    virtual void formDensity() = 0;

  }; // class Quantum

}; // namespace ChronusQ

#endif
