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
#include <cqlinalg/blas1.hpp>

// Debug print (triggers WaveFunction, etc)
//#define _QuantumDebug

namespace ChronusQ {

  enum DENSITY_TYPE {
    SCALAR=0,MZ,MY,MX
  }; ///< Enumerate the types of densities for contraction

  // Helper function for operator traces
  // see src/quantum/properties.cxx for docs

  template <typename RetTyp, typename Left, typename Right>
  static inline RetTyp OperatorTrace(size_t N, const Left& op1 , 
    const Right& op2) {
    return InnerProd<RetTyp>(N,op1,1,op2,1);
  } 

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

    // Helper functions for the automatic evaluation of properties
    // see include/quantum/properties.hpp for documentation

    template <typename RetTyp, DENSITY_TYPE DenTyp, typename Op>
    RetTyp OperatorSpinCombine(const Op&);

  public:

    CQMemManager& memManager; ///< Memory manager for matrix allocation

    int   nC;   ///< Number of spin components
    bool  iCS;  ///< is closed shell?

    // 1PDM storage

    oper_t_coll onePDM;  ///< 1PDM array (Scalar + Magnetization)

    // Property storage

    // Length gauge electric multipoles
    cartvec_t elecDipole;     ///< Electric Dipole in the length gauge
    cartmat_t elecQuadrupole; ///< Electric Quadrupole in the length gauge
    cartrk3_t elecOctupole;   ///< Electric Octupole in the length gauge
    
    // Spin expectation values
    cartvec_t SExpect; ///< Expectation values of Sx, Sy and Sz
    double    SSq;     ///< Expectation value of S^2

    // Energy expectation values
    double OBEnergy;
    double MBEnergy;
    double totalEnergy;

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
      size_t N = 0, bool doAlloc = true): 
      memManager(mem), nC(_nC), iCS(_iCS), 
      elecDipole({0.,0.,0.}),
      elecQuadrupole{{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
      elecOctupole{
        {
          {{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
          {{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
          {{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}}
        }}, SExpect({0.,0.,0.}), SSq(0.), OBEnergy(0.), MBEnergy(0.),
      totalEnergy(0.) {

        // Allocate densities
        if( N != 0 and doAlloc ) alloc(N);

    }; // Quantum::Quantum
    
       
    // See include/quantum/impl.hpp for documentation 
    // on the following constructors

    // Different type
    template <typename U> Quantum(const Quantum<U> &, int dummy = 0);
    template <typename U> Quantum(Quantum<U> &&     , int dummy = 0);

    // Same type
    Quantum(const Quantum &);
    Quantum(Quantum &&);     

    /**
     *  Deconstructor
     */ 
    ~Quantum(){ dealloc(); }

    // Member functions

    // Deallocation (see include/quantum/impl.hpp for docs)
    void alloc(size_t);
    void dealloc();


    // Public interfaces for property evaluation
      
    /**
     *  \brief Computes a 1-body property through a trace with the
     *  proper components of the 1PDM.
     *
     *  \param [template] DenTyp Which spin component of the 1PDM to trace with
     *  \param [in]       op     Square matrix to trace with 1PDM
     *  \returns          Trace of op with the DenTyp 1PDM (cast to type RetTyp)
     */ 
    template <typename RetTyp, DENSITY_TYPE DenTyp, typename Op>
    inline RetTyp computeOBProperty(const Op &op) {
      return OperatorSpinCombine<RetTyp,DenTyp>(op);
    }; // Quantum<T>::computeOBProperty (single operator)

    /**
     *
     */ 
    template <typename RetTyp, DENSITY_TYPE DenTyp, typename Op>
    inline std::vector<RetTyp> computeOBProperty(const std::vector<Op> &opv) {
      std::vector<RetTyp> results;
      for(auto &op : opv) 
        results.emplace_back(computeOBProperty<RetTyp,DenTyp>(op));
      return results;
    }; // Quantum<T>::computeOBProperty (many operators)

    // Procedural
      
    /**
     *  Function to form the density. Pure virtual
     */ 
    virtual void formDensity() = 0;

    /**
     *  Function to compute the energy expectation value(s)
     */ 
    virtual void computeEnergy() = 0;

  }; // class Quantum

}; // namespace ChronusQ

#endif
