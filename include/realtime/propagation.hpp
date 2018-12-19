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
#ifndef __INCLUDED_REALTIME_PROPAGATION_HPP__
#define __INCLUDED_REALTIME_PROPAGATION_HPP__

#include <realtime.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <chronusqefp.hpp>

#include <util/matout.hpp>
#include <unsupported/Eigen/MatrixFunctions>

template <size_t N, typename T>
std::array<T,N> valarray2array(const std::valarray<T> &x) {
  assert( x.size() == N );

  std::array<T,N> arr;
  std::copy_n(&x[0],N,arr.begin());
  return arr;
 
};

namespace ChronusQ {

  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::doPropagation(EFPBase* EFP_,bool EFP_bool) {

    printRTHeader();
    double wf_denp = 0.;
    
    bool Start(false); // Start the MMUT iterations
    bool FinMM(false); // Wrap up the MMUT iterations

    size_t NB = propagator_.aoints.basisSet().nBasis;

    for( curState.xTime = 0., curState.iStep = 0; 
         curState.xTime <= (intScheme.tMax + intScheme.deltaT/4); 
         curState.xTime += intScheme.deltaT, curState.iStep++ ) {

      // Perturbation for the current time
      EMPerturbation pert_t = pert.getPert(curState.xTime);



#if 1
      // Determine the step type for the current integration step 
      if( intScheme.intAlg == MMUT ) {

        // "Start" the MMUT if this is the first step or we just
        // "Finished" the MMUT segment
        Start = ( curState.iStep == 0 ) or FinMM;

        // "Start" the MMUT if the current step index is a restart
        // step
        if( intScheme.iRstrt > 0 ) 
          Start = Start or ( curState.iStep % intScheme.iRstrt == 0 );

        // "Finish" the MMUT if this is the last step
        FinMM = ( curState.iStep == (intScheme.tMax / intScheme.deltaT) );

        // TODO: "Finish" the MMUT if the field turns on or off
        FinMM = FinMM or curState.iStep < 2;
          
        // "Finish" the MMUT if the next step will be a restart step
        if( intScheme.iRstrt > 0 ) 
          FinMM = FinMM or ( (curState.iStep + 1) % intScheme.iRstrt == 0 );

        // If "Starting" or "Finishing" the MMUT, the step type is
        // the specified restart step type, else it is the MMUT step
        if( Start or FinMM ) curState.curStep = intScheme.rstStep;
        else                 curState.curStep = ModifiedMidpoint;

        if( Start or FinMM ) std::cout << "  *** Restarting MMUT ***\n";

      // For non leapfrog scheme, the step type is constant
      } else if ( intScheme.intAlg == ExpMagnus2 )
        curState.curStep = ExplicitMagnus2;
#else
      curState.curStep = ForwardEuler;
#endif












      // Handle density copies / swaps for the current step
      //  + Determine the step size
        
      if( curState.curStep == ModifiedMidpoint ) {
        // Swap the saved density with the SingleSlater density
          
        // DOSav(k) = DO(k)
        // DO(k)    = DO(k-1)
        for(auto i = 0; i < DOSav.size(); i++)
          Swap(memManager_.template getSize<dcomplex>(DOSav[i]),
            DOSav[i],1,propagator_.onePDMOrtho[i],1);

        curState.stepSize = 2. * intScheme.deltaT;

      } else {
        // Save a copy of the SingleSlater density in the saved density
        // storage 
          
        // DOSav(k) = DO(k)
        for(auto i = 0; i < DOSav.size(); i++)
          std::copy_n(propagator_.onePDMOrtho[i],
            memManager_.template getSize<dcomplex>(DOSav[i]),
            DOSav[i]);

     
        curState.stepSize = intScheme.deltaT;

      }





      
      // Form the Fock matrix at the current time
      formFock(false,EFP_,EFP_bool,curState.xTime);

      auto EFP_1 = dynamic_cast< EFP<IntsT,dcomplex>* >(EFP_);
      if(EFP_bool == true and EFP_1 != NULL){
        EFP_1->Wavefunction_dependent(&wf_denp);  
      }

      // Compute properties for D(k) 
      propagator_.computeProperties(pert_t,wf_denp);
      std::cout << "EFP_ind:" << wf_denp << std::endl;

      data.Time.push_back(curState.xTime);
      data.Energy.push_back(propagator_.totalEnergy);
      data.ElecDipole.push_back(propagator_.elecDipole);
      if( pert_t.fields.size() > 0 )
      data.ElecDipoleField.push_back( pert_t.getDipoleAmp(Electric) );


      // Print progress line in the output file
      printRTStep();




      // Orthonormalize the AO Fock matrix
      // F(k) -> FO(k)
      propagator_.ao2orthoFock();


      // Form the propagator from the orthonormal Fock matrix
      // FO(k) -> U**H(k) = exp(- i * dt * FO(k) )
      formPropagator();

      // Propagator the orthonormal density matrix
      // DO (in propagator_) will now store DO(k+1)
      //
      // DO(k+1) = U**H(k) * DO * U(k)
      // - Where DO is what is currently stored in propagator_
      //
      // ***
      // This function also transforms DO(k+1) to the AO
      // basis ( DO(k+1) -> D(k+1) in propagator_ ) and
      // computes the change in density from the previous 
      // AO density ( delD = D(k+1) - D(k) ) 
      // ***
      propagateWFN();

    } // Time loop

  //mathematicaPrint(std::cerr,"Dipole-X",&data.ElecDipole[0][0],
  //  curState.iStep,1,curState.iStep,3);

    if( savFile.exists() ) {
      savFile.safeWriteData("RT/TIME",&data.Time[0],{data.Time.size()});
      savFile.safeWriteData("RT/ENERGY",&data.Energy[0],{data.Time.size()});
      savFile.safeWriteData("RT/LEN_ELEC_DIPOLE",&data.ElecDipole[0][0],
        {data.Time.size(),3});

      if( data.ElecDipoleField.size() > 0 )
      savFile.safeWriteData("RT/LEN_ELEC_DIPOLE_FIELD",&data.ElecDipoleField[0][0],
        {data.Time.size(),3});
    }

  }; // RealTime::doPropagation


  /**
   *  \brief Form the adjoint of the unitary propagator
   *
   *  \f[
   *    U = \exp\left( -i \delta t F \right) 
   *      = \exp\left( -\frac{i\delta t}{2} 
   *                    \left(F^S \otimes I_2 + F^k \sigma_k\right) \right) 
   *      = \frac{1}{2}U^S \otimes I_2 + \frac{1}{2} U^k \otimes \sigma_k
   *  \f]
   */ 
  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::formPropagator() {

    size_t NB = propagator_.aoints.basisSet().nBasis;
    
    // Form U

    // Restricted
    if( UH.size() == 1 ) {
      // See docs for factor of 2
      MatExp('D',NB,dcomplex(0.,-curState.stepSize/2.),
        propagator_.fockMatrixOrtho[SCALAR],NB,UH[SCALAR],NB,memManager_);

      Scale(NB*NB,dcomplex(2.),UH[SCALAR],1);

    // Unrestricted
    } else if( UH.size() == 2 ) {

      // Transform SCALAR / MZ -> ALPHA / BETA
      for(auto i = 0; i < NB*NB; i++) {
        dcomplex tmp = propagator_.fockMatrixOrtho[SCALAR][i];

        propagator_.fockMatrixOrtho[SCALAR][i] = 
          0.5 * (tmp + propagator_.fockMatrixOrtho[MZ][i]);

        propagator_.fockMatrixOrtho[MZ][i] = 
          0.5 * (tmp - propagator_.fockMatrixOrtho[MZ][i]);

      }

      MatExp('D',NB,dcomplex(0.,-curState.stepSize),
        propagator_.fockMatrixOrtho[SCALAR],NB,UH[SCALAR],NB,memManager_);
      MatExp('D',NB,dcomplex(0.,-curState.stepSize),
        propagator_.fockMatrixOrtho[MZ],NB,UH[MZ],NB,memManager_);

      // Transform ALPHA / BETA -> SCALAR / MZ
      for(auto i = 0; i < NB*NB; i++) {
        dcomplex tmp = UH[SCALAR][i];

        UH[SCALAR][i] = tmp + UH[MZ][i];
        UH[MZ][i]     = tmp - UH[MZ][i];

      }

    // Generalized (2C)
    } else {

      dcomplex *SCR  = memManager_.malloc<dcomplex>(8*NB*NB);
      dcomplex *F2C  = SCR;
      dcomplex *UH2C = F2C + 4*NB*NB;

      SpinGather(NB,F2C,2*NB,propagator_.fockMatrixOrtho[SCALAR],NB,
        propagator_.fockMatrixOrtho[MZ],NB,propagator_.fockMatrixOrtho[MY],NB,
        propagator_.fockMatrixOrtho[MX],NB);

      MatExp('D',2*NB,dcomplex(0.,-curState.stepSize),F2C,2*NB,UH2C,2*NB,
        memManager_);

      SpinScatter(NB,UH2C,2*NB,UH[SCALAR],NB,UH[MZ],NB,UH[MY],NB,UH[MX],NB);

      memManager_.free(SCR);
    }

#if 0

    prettyPrintSmart(std::cout,"UH Scalar",UH[SCALAR],NB,NB,NB);
    if( UH.size() > 1 )
      prettyPrintSmart(std::cout,"UH MZ",UH[MZ],NB,NB,NB);

#endif
    
  }; // RealTime::formPropagator



  
  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::propagateWFN() {

    size_t NB = propagator_.aoints.basisSet().nBasis;
    size_t NC = propagator_.nC;
    dcomplex *SCR  = memManager_.template malloc<dcomplex>(NC*NC*NB*NB);
    dcomplex *SCR1 = memManager_.template malloc<dcomplex>(NC*NC*NB*NB);

    if( UH.size() <= 2 ) {

      // Create X(S) = (U**H * DO)(S) in SCR 

      // SCR = 0.5 * U(S)**H * DO(S)
      Gemm('N','N',NB,NB,NB,dcomplex(0.5),UH[SCALAR],NB,
        propagator_.onePDMOrtho[SCALAR],NB,dcomplex(0.),SCR,NB);

      // SCR += 0.5 * U(Z)**H * DO(Z)
      if( UH.size() != 1 )
        Gemm('N','N',NB,NB,NB,dcomplex(0.5),UH[MZ],NB,
          propagator_.onePDMOrtho[MZ],NB,dcomplex(1.),SCR,NB);





      // Create X(Z) = (U**H * DO)(Z) in SCR1
        
      if( propagator_.onePDMOrtho.size() > 1 ) {

        // SCR1 = 0.5 * U(S)**H * DO(Z)
        Gemm('N','N',NB,NB,NB,dcomplex(0.5),UH[SCALAR],NB,
          propagator_.onePDMOrtho[MZ],NB,dcomplex(0.),SCR1,NB);

        // SCR1 += 0.5 * U(Z)**H * DO(S)
        if( UH.size() != 1 )
          Gemm('N','N',NB,NB,NB,dcomplex(0.5),UH[MZ],NB,
            propagator_.onePDMOrtho[SCALAR],NB,dcomplex(1.),SCR1,NB);

      }






      // DO(S) = 0.5 * ( X(S) * U(S) + X(Z) * U(Z) )
      //       = 0.5 * ( SCR  * U(S) + SCR1 * U(Z) )

      // DO(S) = 0.5 * SCR * U(S)
      Gemm('N','C',NB,NB,NB,dcomplex(0.5),SCR,NB,UH[SCALAR],NB,dcomplex(0.),
        propagator_.onePDMOrtho[SCALAR],NB);
 
      // DO(S) += 0.5 * SCR1 * U(Z)
      if( UH.size() != 1 )
        Gemm('N','C',NB,NB,NB,dcomplex(0.5),SCR1,NB,UH[MZ],NB,dcomplex(1.),
          propagator_.onePDMOrtho[SCALAR],NB);


      if( UH.size() != 1) {

        // DO(Z) = 0.5 * ( X(S) * U(Z) + X(Z) * U(S) )
        //       = 0.5 * ( SCR  * U(Z) + SCR1 * U(S) )
          
        // DO(Z) = 0.5 * SCR * U(Z)
        Gemm('N','C',NB,NB,NB,dcomplex(0.5),SCR,NB,UH[MZ],NB,dcomplex(0.),
          propagator_.onePDMOrtho[MZ],NB);
 
        // DO(Z) += 0.5 * SCR1 * U(S)
        Gemm('N','C',NB,NB,NB,dcomplex(0.5),SCR1,NB,UH[SCALAR],NB,dcomplex(1.),
          propagator_.onePDMOrtho[MZ],NB);

      }

    } else {

      // Gather DO
      dcomplex *DO = memManager_.malloc<dcomplex>(4*NB*NB);
      SpinGather(NB,DO,2*NB,propagator_.onePDMOrtho[SCALAR],NB,
        propagator_.onePDMOrtho[MZ],NB,propagator_.onePDMOrtho[MY],NB,
        propagator_.onePDMOrtho[MX],NB);

      // Gather UH into SCR
      SpinGather(NB,SCR,2*NB,UH[SCALAR],NB,UH[MZ],NB,UH[MY],NB,UH[MX],NB);

      // SCR1 = U**H * DO
      Gemm('N','N',2*NB,2*NB,2*NB,dcomplex(1.),SCR,2*NB,DO,2*NB,dcomplex(0.),
        SCR1,2*NB);

      // DO = SCR1 * U
      Gemm('N','C',2*NB,2*NB,2*NB,dcomplex(1.),SCR1,2*NB,SCR,2*NB,dcomplex(0.),
        DO,2*NB);

      // Scatter DO
      SpinScatter(NB,DO,2*NB,propagator_.onePDMOrtho[SCALAR],NB,
        propagator_.onePDMOrtho[MZ],NB,propagator_.onePDMOrtho[MY],NB,
        propagator_.onePDMOrtho[MX],NB);

      // Free memory
      memManager_.free(DO);
    }


    propagator_.ortho2aoDen();

    memManager_.free(SCR,SCR1);

  }; // RealTime::propagatorWFN




}; // namespace ChronusQ

#endif
