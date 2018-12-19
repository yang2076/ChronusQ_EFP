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
#ifndef __INCLUDED_SINGLESLATER_QUANTUM_HPP__
#define __INCLUDED_SINGLESLATER_QUANTUM_HPP__
//#define _DEBUGGIAO 

#include <singleslater.hpp>
#include <chronusqefp.hpp>
#include <cqlinalg/blasext.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/blas3.hpp>
#include <quantum/properties.hpp>

namespace ChronusQ {

  /**
   *  \brief Forms the 1PDM using a set of orbitals 
   *
   *  specialization of Quantum::formDensity. Populates / overwrites
   *  onePDM storage
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::formDensity() {

    size_t NB  = this->aoints.basisSet().nBasis * nC;
    size_t NB2 = NB*NB;


    // Form the 1PDM on the root MPI process as slave processes
    // do not posses the up-to-date MO coefficients
    if( MPIRank(comm) == 0 ) {

      if(nC == 1) { 

        // DS = DA = CA * CA**H
        Gemm('N', 'C', NB, NB, this->nOA, MatsT(1.), this->mo1, NB, this->mo1,
          NB, MatsT(0.), this->onePDM[SCALAR], NB);

        if(not iCS) {

          // DZ = DB = CB * CB**H
          Gemm('N', 'C', NB, NB, this->nOB, MatsT(1.), this->mo2, NB, 
            this->mo2, NB, MatsT(0.), this->onePDM[MZ], NB);

          // DS = DA + DB
          // DZ = DA - DB
          for(auto j = 0; j < NB2; j++) {
            MatsT tmp = this->onePDM[SCALAR][j];

            this->onePDM[SCALAR][j] = 
              this->onePDM[SCALAR][j] + this->onePDM[MZ][j]; 

            this->onePDM[MZ][j]     = tmp - this->onePDM[MZ][j]; 
          }

        } else {

          // DS = 2 * DA
          std::transform(this->onePDM[SCALAR], this->onePDM[SCALAR] + NB2,
            this->onePDM[SCALAR],[](MatsT a){ return 2.*a; }
          );

        }
      } else {

        MatsT * SCR = this->memManager.template malloc<MatsT>(NB2);

        Gemm('N', 'C', NB, NB, this->nO, MatsT(1.), this->mo1, NB, this->mo1,
          NB, MatsT(0.), SCR, NB);

        SpinScatter(NB/2,SCR,NB,this->onePDM[SCALAR],NB/2,this->onePDM[MZ],
          NB/2,this->onePDM[MY],NB/2,this->onePDM[MX],NB/2);

        this->memManager.free(SCR);

      }

    }


#ifdef CQ_ENABLE_MPI

    // Broadcast the 1PDM to all MPI processes
    if( MPISize(comm) > 1 ) {
      std::cerr  << "  *** Scattering the 1PDM ***\n";
      for(auto k = 0; k < this->onePDM.size(); k++)
        MPIBCast(this->onePDM[k],NB*NB/nC/nC,0,comm);
    }

#endif


#if 0
      print1PDM(std::cerr);
#endif

  }; // SingleSlater<T>::formDensity


  /**
   *  \brief Computes the total field free energy of a single slater determinent
   *
   *  Given a 1PDM and a Fock matrix (specifically the core Hamiltonian and
   *  G[D]), compute the energy.
   *
   *  \warning Assumes that density and Fock matrix have the appropriate form
   *
   *  F(S) = F(A) + F(B)
   *  F(Z) = F(A) - F(B)
   *  ...
   *
   *  \f[
   *     E = \frac{1}{2} \mathrm{Tr}[\mathbf{P}(\mathbf{H} + \mathbf{F})] +
   *     V_{NN}
   *  \f]
   *
   *  Specialization of Quantum<T>::computeEnergy, populates / overwrites 
   *  OBEnergy and MBEnergy
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeEnergy() {

    ROOT_ONLY(comm);

    // Scalar core hamiltonian contribution to the energy
    this->OBEnergy = 
      this->template computeOBProperty<double,DENSITY_TYPE::SCALAR>(
         coreH[SCALAR]);

    
    
    // One body Spin Orbit
    dcomplex SOEnergy(0.,0.);
    if(coreH.size() > 1) {
      SOEnergy = 
        this->template computeOBProperty<dcomplex,DENSITY_TYPE::MZ>(
           coreH[MZ]);
      SOEnergy += 
        this->template computeOBProperty<dcomplex,DENSITY_TYPE::MY>(
           coreH[MY]);
      SOEnergy += 
        this->template computeOBProperty<dcomplex,DENSITY_TYPE::MX>(
           coreH[MX]);
    };

 
    this->OBEnergy += std::real(SOEnergy);


    this->OBEnergy *= 0.5;

    // Compute many-body contribution to energy
    // *** These calls are safe as proper zeros are returned by
    // property engine ***
    this->MBEnergy = 
      this->template computeOBProperty<double,DENSITY_TYPE::SCALAR>(this->twoeH[SCALAR]);
    this->MBEnergy += 
      this->template computeOBProperty<double,DENSITY_TYPE::MZ>(this->twoeH[MZ]);
    this->MBEnergy += 
      this->template computeOBProperty<double,DENSITY_TYPE::MY>(this->twoeH[MY]);
    this->MBEnergy += 
      this->template computeOBProperty<double,DENSITY_TYPE::MX>(this->twoeH[MX]);

    this->MBEnergy *= 0.25;
    this->NREnergy = this->aoints.molecule().nucRepEnergy;


    // Assemble total energy
    this->totalEnergy = 
      this->OBEnergy + this->MBEnergy + this->NREnergy + this->EFPEnergy;

    // Sanity checks
    assert( not std::isnan(this->OBEnergy) );
    assert( not std::isnan(this->MBEnergy) );
    assert( not std::isnan(this->EFPEnergy) );
    assert( not std::isinf(this->OBEnergy) );
    assert( not std::isinf(this->MBEnergy) );
    assert( not std::isinf(this->EFPEnergy) );

  }; // SingleSlater<T>::computeEnergy

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeMultipole(EMPerturbation &pert) {

    ROOT_ONLY(comm);
    // Compute elecric contribution to the dipoles
    for(auto iXYZ = 0; iXYZ < 3; iXYZ++) 
      this->elecDipole[iXYZ] = -this->template computeOBProperty<double,SCALAR>(this->aoints.lenElecDipole[iXYZ]);

    // Nuclear contributions to the dipoles
    for(auto &atom : this->aoints.molecule().atoms)
      MatAdd('N','N',3,1,1.,&this->elecDipole[0],3,double(atom.atomicNumber),
        &atom.coord[0],3,&this->elecDipole[0],3);


    // Electric contribution to the quadrupoles
    for(size_t iXYZ = 0, iX = 0; iXYZ < 3; iXYZ++)
    for(size_t jXYZ = iXYZ     ; jXYZ < 3; jXYZ++, iX++){

      this->elecQuadrupole[iXYZ][jXYZ] = -
        this->template computeOBProperty<double,SCALAR>(this->aoints.lenElecQuadrupole[iX]);
      
      this->elecQuadrupole[jXYZ][iXYZ] = this->elecQuadrupole[iXYZ][jXYZ]; 
    }
    
    // Nuclear contributions to the quadrupoles
    for(auto &atom : this->aoints.molecule().atoms)
    for(size_t iXYZ = 0; iXYZ < 3; iXYZ++)
    for(size_t jXYZ = 0; jXYZ < 3; jXYZ++) 
      this->elecQuadrupole[iXYZ][jXYZ] +=
        atom.atomicNumber * atom.coord[iXYZ] * atom.coord[jXYZ];

    

    // Electric contribution to the octupoles
    for(size_t iXYZ = 0, iX = 0; iXYZ < 3; iXYZ++)
    for(size_t jXYZ = iXYZ     ; jXYZ < 3; jXYZ++)
    for(size_t kXYZ = jXYZ     ; kXYZ < 3; kXYZ++, iX++){

      this->elecOctupole[iXYZ][jXYZ][kXYZ] = -
        this->template computeOBProperty<double,SCALAR>(
          this->aoints.lenElecOctupole[iX]);

      this->elecOctupole[iXYZ][kXYZ][jXYZ] = this->elecOctupole[iXYZ][jXYZ][kXYZ]; 

      this->elecOctupole[jXYZ][iXYZ][kXYZ] = this->elecOctupole[iXYZ][jXYZ][kXYZ]; 

      this->elecOctupole[jXYZ][kXYZ][iXYZ] = this->elecOctupole[iXYZ][jXYZ][kXYZ]; 

      this->elecOctupole[kXYZ][iXYZ][jXYZ] = this->elecOctupole[iXYZ][jXYZ][kXYZ]; 

      this->elecOctupole[kXYZ][jXYZ][iXYZ] = this->elecOctupole[iXYZ][jXYZ][kXYZ]; 
    }

    // Nuclear contributions to the octupoles
    for(auto &atom : this->aoints.molecule().atoms)
    for(size_t iXYZ = 0; iXYZ < 3; iXYZ++)
    for(size_t jXYZ = 0; jXYZ < 3; jXYZ++)
    for(size_t kXYZ = 0; kXYZ < 3; kXYZ++)
      this->elecOctupole[iXYZ][jXYZ][kXYZ] +=
        atom.atomicNumber * atom.coord[iXYZ] * atom.coord[jXYZ] *
        atom.coord[kXYZ];

  };


  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeSpin() {

    ROOT_ONLY(comm);

    this->SExpect[0] = 0.5 * this->template computeOBProperty<double,MX>(
      this->aoints.overlap);
    this->SExpect[1] = 0.5 * this->template computeOBProperty<double,MY>(
      this->aoints.overlap);
    this->SExpect[2] = 0.5 * this->template computeOBProperty<double,MZ>(
      this->aoints.overlap);

    if( this->onePDM.size() == 1 ) this->SSq = 0;
    else {
      size_t NB = this->aoints.basisSet().nBasis;
      MatsT * SCR  = this->memManager.template malloc<MatsT>(NB*NB);
      MatsT * SCR2 = this->memManager.template malloc<MatsT>(NB*NB);


      // SCR2 = S * D(S) * S
/*      
      Gemm('N','C',NB,NB,NB,T(1.),aoints.overlap,NB,this->onePDM[SCALAR],NB,
        T(0.),SCR,NB);
      Gemm('N','C',NB,NB,NB,T(1.),SCR,NB,aoints.overlap,NB,T(0.),SCR2,NB);
*/
      Gemm('N','N',NB,NB,NB,MatsT(1.),this->aoints.overlap,NB,this->onePDM[SCALAR],NB,
        MatsT(0.),SCR,NB);
      Gemm('N','C',NB,NB,NB,MatsT(1.),this->aoints.overlap,NB,SCR,NB,MatsT(0.),SCR2,NB);

      
      this->SSq = 3 * this->nO - (3./2.) * 
        this->template computeOBProperty<double,SCALAR>(SCR2);
  

      // SCR2 = D(Z) * S * D(Z)
      Gemm('N','N',NB,NB,NB,MatsT(1.),this->aoints.overlap,NB,this->onePDM[MZ],NB,
        MatsT(0.),SCR,NB);
      Gemm('N','C',NB,NB,NB,MatsT(1.),this->aoints.overlap,NB,SCR,NB,MatsT(0.),SCR2,NB);
      
      this->SSq += 0.5 * this->template computeOBProperty<double,MZ>(SCR2);

      if( this->onePDM.size() > 2 ) {
  
        // SCR2 = D(Y) * S * D(Y)
        Gemm('N','N',NB,NB,NB,MatsT(1.),this->aoints.overlap,NB,this->onePDM[MY],NB,
          MatsT(0.),SCR,NB);
        Gemm('N','C',NB,NB,NB,MatsT(1.),this->aoints.overlap,NB,SCR,NB,MatsT(0.),SCR2,NB);
        
        this->SSq += 0.5 * this->template computeOBProperty<double,MY>(SCR2);


        // SCR2 = D(X) * S * D(X)
        Gemm('N','N',NB,NB,NB,MatsT(1.),this->aoints.overlap,NB,this->onePDM[MX],NB,
          MatsT(0.),SCR,NB);
        Gemm('N','C',NB,NB,NB,MatsT(1.),this->aoints.overlap,NB,SCR,NB,MatsT(0.),SCR2,NB);
        
        this->SSq += 0.5 * this->template computeOBProperty<double,MX>(SCR2);

      }

      for(auto i = 0; i < 3; i++) 
        this->SSq += 4 * this->SExpect[i] * this->SExpect[i];

      this->SSq *= 0.25;

      this->memManager.free(SCR,SCR2);
    }

  };
  

}; // namespace ChronusQ

#endif
