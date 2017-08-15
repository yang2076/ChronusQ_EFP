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

#include <aointegrals.hpp>
#include <cqlinalg.hpp>
#include <cqlinalg/svd.hpp>
#include <cqlinalg/blasutil.hpp>

#include <physcon.hpp>
#include <util/matout.hpp>


// Include inhouse integral builder
#include "aointegrals_builders_inhouse.cxx"

// X2C_DEBUG_LEVEL == 1 - Common problem areas
// X2C_DEBUG_LEVEL == 2 - 1 + Intermediates (W, 4C Core Hamiltonian, etc)
// X2C_DEBUG_LEVEL >= 3 - Print EVERYTHING
#ifndef X2C_DEBUG_LEVEL
#  define X2C_DEBUG_LEVEL 0
#endif

namespace ChronusQ {


  /**
   *  \brief Compute the X2C Core Hamiltonian
   */ 
  void AOIntegrals::computeX2CCH(std::vector<double*> &CH) {

    size_t NP = basisSet_.nPrimitive;
    size_t NB = basisSet_.nBasis;

    // Transformation matrix
    double *UK = memManager_.malloc<double>(NP*NP);

    // Allocate Scratch Space (enough for 2*NP x 2*NP complex matricies)
    double   *SCR1  = memManager_.malloc<double>(8*NP*NP);
    dcomplex *CSCR1 = reinterpret_cast<dcomplex*>(SCR1);


    // Uncontract the basis
    auto uncontractedShells = basisSet_.uncontractShells();

    // Compute S and T integrals in the unonctracted basis
    auto _overlap = OneEDriver(libint2::Operator::overlap,uncontractedShells);
    auto _kinetic = OneEDriver(libint2::Operator::kinetic,uncontractedShells);


    // Compute V + PVP integrals
#if 1
    auto _potential = OneEDriverLocal<1,true>(
                std::bind(
                  static_cast<
                    std::vector<std::vector<double>>
                    (AOIntegrals::*)(
                      libint2::ShellPair&,libint2::Shell&,libint2::Shell&
                    )
                > (&AOIntegrals::computePotentialV),this,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                uncontractedShells);

    auto _SL = OneEDriverLocal<3,false>(
                std::bind(
                  static_cast<
                    std::vector<std::vector<double>>
                    (AOIntegrals::*)(
                      libint2::ShellPair&,libint2::Shell&,libint2::Shell&
                    )
                > (&AOIntegrals::computeSL),this,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                uncontractedShells);

    auto _PVdP = OneEDriverLocal<1,true>(
                std::bind(
                  static_cast<
                    std::vector<std::vector<double>>
                    (AOIntegrals::*)(
                      libint2::ShellPair&,libint2::Shell&,libint2::Shell&
                    )
                > (&AOIntegrals::computepVdotp),this,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                uncontractedShells);
#else
    auto _potential = OneEDriverLocal<1,true>(
                std::bind(
                  static_cast<
                    std::vector<std::vector<double>>
                    (AOIntegrals::*)(
                      const shell_set &,
                      libint2::ShellPair&,libint2::Shell&,libint2::Shell&
                    )
                > (&AOIntegrals::computePotentialV),this,molecule_.chargeDist,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                uncontractedShells);

    auto _SL = OneEDriverLocal<3,false>(
                std::bind(
                  static_cast<
                    std::vector<std::vector<double>>
                    (AOIntegrals::*)(
                      const shell_set &,
                      libint2::ShellPair&,libint2::Shell&,libint2::Shell&
                    )
                > (&AOIntegrals::computeSL),this,molecule_.chargeDist,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                uncontractedShells);

    auto _PVdP = OneEDriverLocal<1,true>(
                std::bind(
                  static_cast<
                    std::vector<std::vector<double>>
                    (AOIntegrals::*)(
                      const shell_set &,
                      libint2::ShellPair&,libint2::Shell&,libint2::Shell&
                    )
                > (&AOIntegrals::computepVdotp),this,molecule_.chargeDist,
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3),
                uncontractedShells);

#endif


    // del -> p
    for(auto &SL : _SL) Scale(NP*NP,-1.,SL,1);

    // Compute the mappings from primitives to CGTOs
    double * mapPrim2Cont = memManager_.malloc<double>(NP*NB);
    basisSet_.makeMapPrim2Cont(_overlap[0],mapPrim2Cont,memManager_);

#if X2C_DEBUG_LEVEL >= 3
    prettyPrintSmart(std::cout,"Overlap",_overlap[0],NP,NP,NP);
    prettyPrintSmart(std::cout,"Kinetic",_kinetic[0],NP,NP,NP);
    prettyPrintSmart(std::cout,"Potential",_potential[0],NP,NP,NP);

    prettyPrintSmart(std::cout,"PV.P",_PVdP[0],NP,NP,NP);
    prettyPrintSmart(std::cout,"(PVxP) X",_SL[0],NP,NP,NP);
    prettyPrintSmart(std::cout,"(PVxP) Y",_SL[1],NP,NP,NP);
    prettyPrintSmart(std::cout,"(PVxP) Z",_SL[2],NP,NP,NP);
#endif

    // Make a copy of the overlap for later
    double * SCPY = memManager_.malloc<double>(NP*NP);
    memcpy(SCPY,_overlap[0],NP*NP*sizeof(double));
    


    // Singular value storage (initially S then T)
    double * SS   = memManager_.malloc<double>(NP);
    
    // Get SVD of uncontracted overlap
    // Store the left singular vectors in S
    SVD('O','N',NP,NP,_overlap[0],NP,SS,reinterpret_cast<double*>(NULL),NP,
      reinterpret_cast<double*>(NULL),NP,memManager_);

    double minSS = *std::min_element(SS,SS+NP);

    if( minSS < 1e-10 ) CErr("Uncontracted Overlap is Singular");

#if X2C_DEBUG_LEVEL >= 1
    prettyPrintSmart(std::cout,"Uncontracted S Singular Values",SS,NP,1,NP);
#endif

    // Form orthonormal transformation matrix in S
    for(auto i = 0ul; i < NP; i++)
      Scale(NP,1./std::sqrt(SS[i]),_overlap[0] + i*NP,1);


    // Transform T into the orthonormal basis
    // T -> TO
    Gemm('T','N',NP,NP,NP,1.,_overlap[0],NP,_kinetic[0],NP,0.,SCR1,NP);
    Gemm('N','N',NP,NP,NP,1.,SCR1,NP,_overlap[0],NP,0.,_kinetic[0],NP);

    // Get the SVD of TO
    // Store the left singular vectors in TO
    SVD('O','N',NP,NP,_kinetic[0],NP,SS,reinterpret_cast<double*>(NULL),NP,
      reinterpret_cast<double*>(NULL),NP,memManager_);

    minSS = *std::min_element(SS,SS+NP);
    if( minSS < 1e-10 )
      CErr("Uncontracted Kinetic Energy Tensor is Singular");

#if X2C_DEBUG_LEVEL >= 3
    prettyPrintSmart(std::cout,"S Transformation",_overlap[0],NP,NP,NP);
    prettyPrintSmart(std::cout,"T Transformation",_kinetic[0],NP,NP,NP);
#endif

    // Form UK = S * T
    Gemm('N','N',NP,NP,NP,1.,_overlap[0],NP,_kinetic[0],NP,0.,UK,NP);

#if X2C_DEBUG_LEVEL >= 2
    prettyPrintSmart(std::cout,"UK",UK,NP,NP,NP);
#endif
    
    // Free up S and T storage as they're no longer used after UK is formed
    memManager_.free(_overlap[0],_kinetic[0]);



    // Allocate and for "P^2" potential
    double *P2P = memManager_.malloc<double>(NP*NP);

    // P2P = UK**T * V * UK
    Gemm('T','N',NP,NP,NP,1.,UK,NP,_potential[0],NP,0.,SCR1,NP);
    Gemm('N','N',NP,NP,NP,1.,SCR1,NP,UK,NP,0.,P2P,NP);

#if X2C_DEBUG_LEVEL >= 3
    prettyPrintSmart(std::cout,"P2P",P2P,NP,NP,NP);
#endif

    // Free up V storage as it's no longer used after P2P is formed 
    memManager_.free(_potential[0]);


    // Transform PVP into the "P^2" basis
    Gemm('T','N',NP,NP,NP,1.,UK,NP,_PVdP[0],NP,0.,SCR1,NP);
    Gemm('N','N',NP,NP,NP,1.,SCR1,NP,UK,NP,0.,_PVdP[0],NP);
    
    // Loop over PVxP terms
    for(auto & SL : _SL ){
      Gemm('T','N',NP,NP,NP,1.,UK,NP,SL,NP,0.,SCR1,NP);
      Gemm('N','N',NP,NP,NP,1.,SCR1,NP,UK,NP,0.,SL,NP);
    }

#if X2C_DEBUG_LEVEL >= 3
    prettyPrintSmart(std::cout,"(P^2) PV.P",_PVdP[0],NP,NP,NP);
    prettyPrintSmart(std::cout,"(P^2) (PVxP) X",_SL[0],NP,NP,NP);
    prettyPrintSmart(std::cout,"(P^2) (PVxP) Y",_SL[1],NP,NP,NP);
    prettyPrintSmart(std::cout,"(P^2) (PVxP) Z",_SL[2],NP,NP,NP);
#endif

    // P^2 -> P^-1
    for(auto i = 0; i < NP; i++) {
      SS[i] = 1./std::sqrt(2*SS[i]);
    }

    // Transform PVP into "P^-1" basis
    for(auto j = 0; j < NP; j++) 
    for(auto i = 0; i < NP; i++){
      _PVdP[0][i + j*NP] *= SS[i] * SS[j];
      for(auto &SL : _SL)
        SL[i + j*NP] *= SS[i] * SS[j];
    } 

#if X2C_DEBUG_LEVEL >= 3
    prettyPrintSmart(std::cout,"(P^-1) PV.P",_PVdP[0],NP,NP,NP);
    prettyPrintSmart(std::cout,"(P^-1) (PVxP) X",_SL[0],NP,NP,NP);
    prettyPrintSmart(std::cout,"(P^-1) (PVxP) Y",_SL[1],NP,NP,NP);
    prettyPrintSmart(std::cout,"(P^-1) (PVxP) Z",_SL[2],NP,NP,NP);
#endif

    // Allocate 4C CORE Hamiltonian

    // CH = [ V    cp       ]
    //      [ cp   W - 2mc^2]
    dcomplex *CH4C = memManager_.malloc<dcomplex>(16*NP*NP);
    memset(CH4C,0,16*NP*NP*sizeof(dcomplex));


    // Allocate W separately  as it's needed later
    size_t LDW = 2*NP;
    dcomplex *W  = memManager_.malloc<dcomplex>(LDW*LDW);
      
    // W = [ W1  W2 ]
    //     [ W3  W4 ]
    dcomplex *W1 = W;
    dcomplex *W2 = W1 + LDW*NP;
    dcomplex *W3 = W1 + NP;
    dcomplex *W4 = W2 + NP;
      
    // W1 = pV.p + i (pVxp)(Z)
    SetMatRE('N',NP,NP,1.,_PVdP[0],NP,W1,LDW);
    SetMatIM('N',NP,NP,1.,_SL[2],  NP,W1,LDW);

    // W4 = conj(W1)
    SetMatRE('N',NP,NP,1., _PVdP[0],NP,W4,LDW);
    SetMatIM('N',NP,NP,-1.,_SL[2],  NP,W4,LDW);
  
    // W2 = (pVxp)(Y) + i (pVxp)(X)
    SetMatRE('N',NP,NP,1.,_SL[1],NP,W2,LDW);
    SetMatIM('N',NP,NP,1.,_SL[0],NP,W2,LDW);

    // W3 = -conj(W2)
    SetMatRE('N',NP,NP,-1.,_SL[1],NP,W3,LDW);
    SetMatIM('N',NP,NP,1., _SL[0],NP,W3,LDW);

#if X2C_DEBUG_LEVEL >= 2
    prettyPrintSmart(std::cout,"W",W,2*NP,2*NP,LDW);
#endif

    // Subtract out 2mc^2 from W diagonals
    double WFact = 2. * SpeedOfLight * SpeedOfLight;
    for(auto j = 0ul; j < 2*NP; j++) W[j + LDW*j] -= WFact;
  


    // Copy W into the 4C CH storage
    dcomplex *CHW = CH4C + 8*NP*NP + 2*NP;
    SetMat('N',2*NP,2*NP,dcomplex(1.),W,LDW,CHW,4*NP);
  

    // Free up PVP memory
    memManager_.free(_PVdP[0]);
    for(auto &SL : _SL ) memManager_.free(SL);


    // P^-1 -> P
    for(auto i = 0; i < NP; i++) SS[i] = 1./SS[i];

    // V = [ P2P  0   ]
    //     [ 0    P2P ]
    dcomplex * V1 = CH4C;
    dcomplex * V2 = V1 + 4*NP*NP + NP;

    SetMatRE('N',NP,NP,1.,P2P,NP,V1,4*NP);
    SetMatRE('N',NP,NP,1.,P2P,NP,V2,4*NP);
    

    // Set the diagonal cp blocks of CH
    // CP = [cp 0  ]
    //      [0  cp ]
    dcomplex *CP11 = CH4C + 8*NP*NP;
    dcomplex *CP12 = CP11 + 4*NP*NP + NP;
    dcomplex *CP21 = CH4C + 2*NP;
    dcomplex *CP22 = CP21 + 4*NP*NP + NP;
   
    for(auto j = 0; j < NP; j++) {
      CP11[j + 4*NP*j] = SpeedOfLight * SS[j];
      CP12[j + 4*NP*j] = SpeedOfLight * SS[j];
      CP21[j + 4*NP*j] = SpeedOfLight * SS[j];
      CP22[j + 4*NP*j] = SpeedOfLight * SS[j];
    }

#if X2C_DEBUG_LEVEL >= 2
    prettyPrintSmart(std::cout,"4C Core Hamiltonian",CH4C,4*NP,4*NP,4*NP);
#endif

    // Diagonalize the 4C CH
    double *CHEV = memManager_.malloc<double>(4*NP);

    HermetianEigen('V','U',4*NP,CH4C,4*NP,CHEV,memManager_);

#if X2C_DEBUG_LEVEL >= 3
    prettyPrintSmart(std::cout,"EV",CHEV,4*NP,1,4*NP);
    prettyPrintSmart(std::cout,"VC",CH4C,4*NP,4*NP,4*NP);
#endif


    // Get pointers to "L" and "S" components of eigenvectors
    dcomplex *L = CH4C + 8*NP*NP;
    dcomplex *S = L + 2*NP;


    // Invert "L"; L -> L^-1
    LUInv(2*NP,L,4*NP,memManager_);




    // Reuse the charge conjugated space for X and Y
    dcomplex *X = CH4C;
    dcomplex *Y = X + 2*NP;

    // Form X = S * L^-1
    Gemm('N','N',2*NP,2*NP,2*NP,dcomplex(1.),S,4*NP,L,4*NP,
      dcomplex(0.),X,4*NP);

    // Form Y = sqrt(1 + X**H * H)
      
    // Y = X**H * H
    Gemm('C','N',2*NP,2*NP,2*NP,dcomplex(1.),X,4*NP,X,4*NP,
      dcomplex(0.),Y,4*NP);

    // Y = Y + I
    for(auto j = 0; j < 2*NP; j++) Y[j + 4*NP*j] += 1.0;

    // Y -> V * y * V**H 
    // XXX: Store the eigenvalues of Y in CHEV
    HermetianEigen('V','U',2*NP,Y,4*NP,CHEV,memManager_);

    // SCR1 -> V * y^-0.25
    for(auto j = 0ul; j < 2*NP; j++)
    for(auto i = 0ul; i < 2*NP; i++)
      CSCR1[i + 2*NP*j] = Y[i + 4*NP*j] * std::pow(CHEV[j],-0.25);

    // Y = SCR1 * SCR1**H
    Gemm('N','C',2*NP,2*NP,2*NP,dcomplex(1.),CSCR1,2*NP,CSCR1,2*NP,
      dcomplex(0.),Y,4*NP);
    
    
#if X2C_DEBUG_LEVEL >= 1
    prettyPrintSmart(std::cout,"X",X,2*NP,2*NP,4*NP);
    prettyPrintSmart(std::cout,"Y",Y,2*NP,2*NP,4*NP);
#endif

    // Build the effective two component CH in "L"
    dcomplex *FullCH2C = L;
 
    // Zero it out
    for(auto j = 0; j < 2*NP; j++)
    for(auto i = 0; i < 2*NP; i++)
      FullCH2C[i + 4*NP*j] = 0.;

    // Copy P2P into spin diagonal blocks of 2C CH
    dcomplex *CH2C1 = FullCH2C;
    dcomplex *CH2C2 = CH2C1 + 4*NP*NP + NP;
    SetMatRE('N',NP,NP,1.,P2P,NP,CH2C1,4*NP);
    SetMatRE('N',NP,NP,1.,P2P,NP,CH2C2,4*NP);

    // Construct 2C CH in the uncontracted basis
    // 2C CH = Y * (V' + cp * X + X**H * cp + X**H * W' * X) * Y

    // SCR1 = cp * X
    for(auto j = 0; j < 2*NP; j++)
    for(auto i = 0; i < NP; i++) {
      CSCR1[i + 2*NP*j] = SpeedOfLight * SS[i] * X[i + 4*NP*j];
      CSCR1[i + NP + 2*NP*j] = SpeedOfLight * SS[i] * X[i + NP + 4*NP*j];
    }

    // 2C CH += SCR1 + SCR1**H
    MatAdd('N','N',2*NP,2*NP,dcomplex(1.),FullCH2C,4*NP,dcomplex(1.),
      CSCR1,2*NP, FullCH2C,4*NP);
    MatAdd('N','C',2*NP,2*NP,dcomplex(1.),FullCH2C,4*NP,dcomplex(1.),
      CSCR1,2*NP, FullCH2C,4*NP);


    // SCR1 = X**H * W
    Gemm('C','N',2*NP,2*NP,2*NP,dcomplex(1.),X,4*NP,W,LDW,dcomplex(0.),
      CSCR1,2*NP);
    
    // 2C CH += SCR1 * X
    Gemm('N','N',2*NP,2*NP,2*NP,dcomplex(1.),CSCR1,2*NP,X,4*NP,
      dcomplex(1.),FullCH2C,4*NP);

    // SCR1 = CH2C * Y
    Gemm('C','N',2*NP,2*NP,2*NP,dcomplex(1.),FullCH2C,4*NP,Y,4*NP,dcomplex(0.),
      CSCR1,2*NP);


    // 2C CH = Y * SCR1
    Gemm('N','N',2*NP,2*NP,2*NP,dcomplex(1.),Y,4*NP,CSCR1,2*NP,
      dcomplex(0.),FullCH2C,4*NP);

#if X2C_DEBUG_LEVEL >= 3
    prettyPrintSmart(std::cout,"2C Core Hamiltonian (P-Space)",
      FullCH2C,2*NP,2*NP,4*NP);
#endif


    // Allocate memory for the uncontracted spin components 
    // of the 2C CH
    dcomplex *HUnS = memManager_.malloc<dcomplex>(NP*NP);
    dcomplex *HUnZ = memManager_.malloc<dcomplex>(NP*NP);
    dcomplex *HUnX = memManager_.malloc<dcomplex>(NP*NP);
    dcomplex *HUnY = memManager_.malloc<dcomplex>(NP*NP);

    SpinScatter(NP,FullCH2C,4*NP,HUnS,NP,HUnZ,NP,HUnY,NP,HUnX,NP);

#if 0
    // XXX: For debugging with old code
    Scale(NP*NP,dcomplex(0.5),HUnS,1);
    Scale(NP*NP,dcomplex(0.5),HUnZ,1);
    Scale(NP*NP,dcomplex(0.5),HUnY,1);
    Scale(NP*NP,dcomplex(0.5),HUnX,1);
#endif

#if X2C_DEBUG_LEVEL >= 3
    prettyPrintSmart(std::cout,"(P) 2C H(S)",HUnS,NP,NP,NP);
    prettyPrintSmart(std::cout,"(P) 2C H(Z)",HUnZ,NP,NP,NP);
    prettyPrintSmart(std::cout,"(P) 2C H(Y)",HUnY,NP,NP,NP);
    prettyPrintSmart(std::cout,"(P) 2C H(X)",HUnX,NP,NP,NP);
#endif
    
    
    // Partition the scratch space into one complex and one real NP x NP 
    // matrix
    double   * SUK   = SCR1;
    dcomplex * CSCR2 = reinterpret_cast<dcomplex*>(SUK + NP*NP);

    // Store the Product of S and UK
    Gemm('N','N',NP,NP,NP,1.,SCPY,NP,UK,NP,0.,SCR1,NP); 

    // Transform the spin components of the 2C CH into R-space
    //
    // H(k) -> SUK * H(k) * (SUK)**H
    //
    // ** Using the fact that H(k) is hermetian
    // CSCR2 = SUK * H(k) -> CSCR2**H = H(k) * (SUK)**H
    // H(k) -> SUK * CSCR2**H
    //

    // Transform H(S)
    Gemm('N','N',NP,NP,NP,dcomplex(1.),SUK,NP,HUnS,NP,dcomplex(0.),
      CSCR2,NP);
    Gemm('N','C',NP,NP,NP,dcomplex(1.),SUK,NP,CSCR2,NP,dcomplex(0.),
      HUnS,NP);

    // Transform H(Z)
    Gemm('N','N',NP,NP,NP,dcomplex(1.),SUK,NP,HUnZ,NP,dcomplex(0.),
      CSCR2,NP);
    Gemm('N','C',NP,NP,NP,dcomplex(1.),SUK,NP,CSCR2,NP,dcomplex(0.),
      HUnZ,NP);

    // Transform H(Y)
    Gemm('N','N',NP,NP,NP,dcomplex(1.),SUK,NP,HUnY,NP,dcomplex(0.),
      CSCR2,NP);
    Gemm('N','C',NP,NP,NP,dcomplex(1.),SUK,NP,CSCR2,NP,dcomplex(0.),
      HUnY,NP);

    // Transform H(X)
    Gemm('N','N',NP,NP,NP,dcomplex(1.),SUK,NP,HUnX,NP,dcomplex(0.),
      CSCR2,NP);
    Gemm('N','C',NP,NP,NP,dcomplex(1.),SUK,NP,CSCR2,NP,dcomplex(0.),
      HUnX,NP);


#if X2C_DEBUG_LEVEL >= 3
    prettyPrintSmart(std::cout,"(R) 2C H(S)",HUnS,NP,NP,NP);
    prettyPrintSmart(std::cout,"(R) 2C H(Z)",HUnZ,NP,NP,NP);
    prettyPrintSmart(std::cout,"(R) 2C H(Y)",HUnY,NP,NP,NP);
    prettyPrintSmart(std::cout,"(R) 2C H(X)",HUnX,NP,NP,NP);
#endif

    // Transform H(k) into the contracted basis

    Gemm('N','N',NB,NP,NP,dcomplex(1.),mapPrim2Cont,NB,HUnS,
      NP,dcomplex(0.),CSCR1,NB);
    Gemm('N','C',NB,NB,NP,dcomplex(1.),mapPrim2Cont,NB,CSCR1,
      NB,dcomplex(0.),HUnS,NB);

    Gemm('N','N',NB,NP,NP,dcomplex(1.),mapPrim2Cont,NB,HUnZ,
      NP,dcomplex(0.),CSCR1,NB);
    Gemm('N','C',NB,NB,NP,dcomplex(1.),mapPrim2Cont,NB,CSCR1,
      NB,dcomplex(0.),HUnZ,NB);

    Gemm('N','N',NB,NP,NP,dcomplex(1.),mapPrim2Cont,NB,HUnY,
      NP,dcomplex(0.),CSCR1,NB);
    Gemm('N','C',NB,NB,NP,dcomplex(1.),mapPrim2Cont,NB,CSCR1,
      NB,dcomplex(0.),HUnY,NB);

    Gemm('N','N',NB,NP,NP,dcomplex(1.),mapPrim2Cont,NB,HUnX,
      NP,dcomplex(0.),CSCR1,NB);
    Gemm('N','C',NB,NB,NP,dcomplex(1.),mapPrim2Cont,NB,CSCR1,
      NB,dcomplex(0.),HUnX,NB);


#if X2C_DEBUG_LEVEL >= 3
    prettyPrintSmart(std::cout,"X2C H(S) (No Scaling)",HUnS,NB,NB,NB);
    prettyPrintSmart(std::cout,"X2C H(Z) (No Scaling)",HUnZ,NB,NB,NB);
    prettyPrintSmart(std::cout,"X2C H(Y) (No Scaling)",HUnY,NB,NB,NB);
    prettyPrintSmart(std::cout,"X2C H(X) (No Scaling)",HUnX,NB,NB,NB);
#endif

    size_t n1, n2;
    for(auto s1(0ul), i(0ul); s1 < basisSet_.nShell; s1++, i+=n1) {
      n1 = basisSet_.shells[s1].size();
      
      size_t L1 = basisSet_.shells[s1].contr[0].l;
      if ( L1 == 0 ) continue;

      size_t Z1 =
        molecule_.atoms[basisSet_.mapSh2Cen[s1]].atomicNumber;


    for(auto s2(0ul), j(0ul); s2 < basisSet_.nShell; s2++, j+=n2) {
      n2 = basisSet_.shells[s2].size();
      
      size_t L2 = basisSet_.shells[s2].contr[0].l;
      if ( L2 == 0 ) continue;

      size_t Z2 =
        molecule_.atoms[basisSet_.mapSh2Cen[s2]].atomicNumber;
    
      dcomplex fudgeFactor = std::sqrt(
        L1 * (L1+1) * (2*L1+1) *
        L2 * (L2+1) * (2*L2+1) / 9. /
        Z1 / Z2
      );

      MatAdd('N','N',n1,n2,dcomplex(1.),HUnZ + i + j*NB,NB,
        -fudgeFactor,HUnZ + i + j*NB,NB, HUnZ + i + j*NB,NB);

      MatAdd('N','N',n1,n2,dcomplex(1.),HUnY + i + j*NB,NB,
        -fudgeFactor,HUnY + i + j*NB,NB, HUnY + i + j*NB,NB);

      MatAdd('N','N',n1,n2,dcomplex(1.),HUnX + i + j*NB,NB,
        -fudgeFactor,HUnX + i + j*NB,NB, HUnX + i + j*NB,NB);

    } // loop s2
    } // loop s1


    GetMatRE('N',NB,NB,1.,HUnS,NB,CH[0],NB);
    GetMatIM('N',NB,NB,1.,HUnZ,NB,CH[1],NB);
    GetMatIM('N',NB,NB,1.,HUnY,NB,CH[2],NB);
    GetMatIM('N',NB,NB,1.,HUnX,NB,CH[3],NB);

#if X2C_DEBUG_LEVEL >= 1
    prettyPrintSmart(std::cout,"Re[ X2C H(S) ]",CH[0],NB,NB,NB);
    prettyPrintSmart(std::cout,"Im[ X2C H(Z) ]",CH[1],NB,NB,NB);
    prettyPrintSmart(std::cout,"Im[ X2C H(Y) ]",CH[2],NB,NB,NB);
    prettyPrintSmart(std::cout,"Im[ X2C H(X) ]",CH[3],NB,NB,NB);
#endif

    memManager_.free(
      UK,
      SCR1 // Scratch space
    );
  };



  /**
   *  \brief Compute the 4C Dirac Hamiltonian
   */ 
  void AOIntegrals::compute4CCH(std::vector<libint2::Shell> &shells,
    double* ) {

  };


}; // namespace ChronusQ
