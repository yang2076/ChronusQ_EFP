/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2018 Li Research Group (University of Washington)
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
#ifndef __INCLUDED_QUANTUM_BASE_HPP__
#define __INCLUDED_QUANTUM_BASE_HPP__

#include <chronusq_sys.hpp>

#include <util/typedefs.hpp>
#include <util/mpi.hpp>

#include <memmanager.hpp>
#include <cqlinalg/blas1.hpp>

#include <fields.hpp>

namespace ChronusQ {

  enum DENSITY_TYPE {
    SCALAR=0,MZ=1,MY=2,MX=3
  }; ///< Enumerate the types of densities for contraction

  // Helper function for operator traces
  // see src/quantum/properties.cxx for docs

  template <typename RetTyp, typename Left, typename Right>
  static inline RetTyp OperatorTrace(size_t N, const Left& op1 , 
    const Right& op2) {
        
    return InnerProd<RetTyp>(N,op1,1,op2,1);

  } 

  /**
   *  \brief The QuantumBase class. The abstraction of information
   *  relating to the Quantum class which are independent of storage
   *  type.
   *
   *  See Quantum for further docs
   *
   */ 
  class QuantumBase {

  protected:
  private:
  public:

    MPI_Comm      comm; ///< MPI Communicator

    CQMemManager& memManager; ///< Memory manager for matrix allocation

    int   nC;   ///< Number of spin components
    bool  iCS;  ///< is closed shell?


    // Property storage

    // Length gauge electric multipoles
    cart_t elecDipole;        ///< Electric Dipole in the length gauge
    cartmat_t elecQuadrupole; ///< Electric Quadrupole in the length gauge
    cartrk3_t elecOctupole;   ///< Electric Octupole in the length gauge
    
    // Spin expectation values
    cart_t SExpect; ///< Expectation values of Sx, Sy and Sz
    double    SSq;  ///< Expectation value of S^2

    // Energy expectation values
    double OBEnergy;   ///< 1-Body operator contribution to the energy
    double MBEnergy;   ///< Many(2)-Body operator contribution to the energy
    double NREnergy;   ///< Nuclear Repulsion Energy
    double EFPEnergy;   ///< EFP impact Energy
    double totalEnergy;///< The total energy



    // Constructors
      
    // Disable default constructor
    QuantumBase() = delete;

    // Default the copy and move constructors
    QuantumBase(const QuantumBase&) = default;
    QuantumBase(QuantumBase&&)      = default;

    /**
     *  QuantumBase Constructor. Constructs a QuantumBase object.
     *
     *  \param [in] mem   CQMemManager to handle to allocation of densities
     *  \param [in] _nC   Number of spin components (1 and 2 are supported)
     *  \param [in] _iCS  Whether or not system is closed shell
     *                    (only used when _nC == 1)
     */ 
    QuantumBase(MPI_Comm c, CQMemManager &mem, size_t _nC, bool _iCS): 
      memManager(mem), nC(_nC), iCS(_iCS), comm(c),
      elecDipole({0.,0.,0.}),
      elecQuadrupole{{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
      elecOctupole{
        {
          {{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
          {{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
          {{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}}
        }}, SExpect({0.,0.,0.}), SSq(0.), OBEnergy(0.), MBEnergy(0.),
      totalEnergy(0.) { }; // Quantum::Quantum



    // Procedural
      
    /**
     *  Function to form the density. Pure virtual
     */ 
    virtual void formDensity() = 0;

    /**
     *  Function to compute the field free energy expectation value(s)
     */ 
    virtual void computeEnergy() = 0;


    /**
     *  Function to compute the energy expectation values including
     *  field terms
     */ 
    void computeEnergy(EMPerturbation &pert, double wf_denp){

      ROOT_ONLY(comm);

      computeEnergy();


      // Field contributions to the energy
      double field_delta(0.);

      // Purely magnetic field contributions here
      auto magDipoleField = pert.getDipoleAmp(Magnetic);


      // Purely electric field contributions here
      auto elecDipoleField = pert.getDipoleAmp(Electric);
      field_delta += 
        elecDipoleField[0] * elecDipole[0] +
        elecDipoleField[1] * elecDipole[1] +
        elecDipoleField[2] * elecDipole[2];

      totalEnergy += field_delta+wf_denp; // Increment total energy

    };

   
    virtual void computeMultipole(EMPerturbation &) = 0;
    virtual void computeSpin() = 0;
    virtual void methodSpecificProperties() = 0;



    inline void computeProperties(EMPerturbation &pert, double wf_denp) {

      ROOT_ONLY(comm);
      
      computeMultipole(pert);
      computeEnergy(pert,wf_denp);
      computeSpin();
      methodSpecificProperties();

    };
    

    // Print functions
    virtual void print1PDM(std::ostream&) = 0;
    void printMultipoles(std::ostream&);
    void printSpin(std::ostream&);
    virtual void printMiscProperties(std::ostream&) = 0;
  }; // class QuantumBase

}; // namespace ChronusQ

#endif
