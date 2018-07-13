#include <aointegrals.hpp>

namespace ChronusQ {

  /**
   *  \brief Compute the ERI of two shell pairs
   *
   *
   *  \param [in] pair1  bra shell pair data for shell1,shell2
   *  \param [in] pair2  ket shell pair data for shell3,shell4
   *  \param [in] shell1
   *  \param [in] shell2
   *  \param [in] shell3
   *  \param [in] shell4
   *
   *  \return ERI of two shell pairs ( shell1 shell2 | shell3 shell4 )
   */
   
  std::vector<double> RealGTOIntEngine::computeERIabcd(libint2::ShellPair &pair1 ,
    libint2::ShellPair &pair2, libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4)  {
    
    double tmpVal=0.0,sqrPQ,PQ;
    std::vector<double> ERI_cart;
    
    int lA[3],lB[3],lC[3],lD[3];

    // compute total angular momentum
    auto lTotal = shell1.contr[0].l + shell2.contr[0].l
                + shell3.contr[0].l + shell4.contr[0].l;

/*
    std::cerr<<"LA "<<shell1.contr[0].l
    <<" LB "<<shell2.contr[0].l
    <<" LC "<<shell3.contr[0].l 
    <<" LD "<<shell4.contr[0].l<<std::endl;
*/

    // pre calculate all the Boys functions 
    // dimension is FmT_2e[shellpair1.prim][shellpair2.prim][lTotal+1]
    std::vector<std::vector<std::vector<double>>> FmT_2e;
    FmT_2e.resize(pair1.primpairs.size());

    double *FmT = new double[lTotal+1];
    int shellpair1_i=0, shellpair2_j ;
    for ( auto &pripair1 : pair1.primpairs ) {
      FmT_2e[shellpair1_i].resize( pair2.primpairs.size() );

      shellpair2_j = 0 ; 
      for ( auto &pripair2 : pair2.primpairs ) {
        sqrPQ = 0.0;
        for ( int mu=0 ; mu<3 ; mu++ ) {
          PQ = ( pripair1.P[mu]-pripair2.P[mu] ); 
          sqrPQ += PQ*PQ;
        }
        auto Zeta = 1.0/pripair1.one_over_gamma;
        auto Eta  = 1.0/pripair2.one_over_gamma;
        
        auto rho = Zeta*Eta/(Zeta+Eta);
        auto T = rho*sqrPQ;
        // calculate Fm(T) list
        computeFmTTaylor( FmT, T, lTotal, 0 );

        for ( int lcurr = 0 ; lcurr < lTotal+1 ; lcurr++ ) {
           if ( std::abs(FmT[lcurr]) < 1.0e-15 ) 
             FmT_2e[shellpair1_i][shellpair2_j].push_back(0.0);
           else
             FmT_2e[shellpair1_i][shellpair2_j].push_back(FmT[lcurr]);
        } // for lcurr
        shellpair2_j ++;
      } // for pripair2
    shellpair1_i++;
    } // for pripair1
    delete[] FmT;

    for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size(); i++) 
    for(int j = 0; j < cart_ang_list[shell2.contr[0].l].size(); j++)
    for(int k = 0; k < cart_ang_list[shell3.contr[0].l].size(); k++)
    for(int l = 0; l < cart_ang_list[shell4.contr[0].l].size(); l++) {
      for (int mu=0 ; mu<3 ; mu++) {
        lA[mu] = cart_ang_list[shell1.contr[0].l][i][mu];
        lB[mu] = cart_ang_list[shell2.contr[0].l][j][mu];
        lC[mu] = cart_ang_list[shell3.contr[0].l][k][mu];
        lD[mu] = cart_ang_list[shell4.contr[0].l][l][mu];
      }  // for mu

      
      tmpVal = twoehRRabcd(pair1,pair2,shell1,shell2,shell3,shell4,
                 FmT_2e,shell1.contr[0].l, lA, shell2.contr[0].l, lB, 
                 shell3.contr[0].l, lC, shell4.contr[0].l, lD ); 
      
      ERI_cart.push_back(tmpVal);
    }   // for l


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) and 
         ( not shell3.contr[0].pure ) and ( not shell4.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return ERI_cart;
    }

    std::vector<double> ERI_sph;

    ERI_sph.assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)
                   *(2*shell3.contr[0].l+1)*(2*shell4.contr[0].l+1)),0.0); 

    cart2sph_2e_transform( shell1.contr[0].l,shell2.contr[0].l,
      shell3.contr[0].l,shell4.contr[0].l,ERI_sph,ERI_cart );

    return ERI_sph; 

  }   // computeERIabcd


  /**
   *  \brief horizontal recursion of ERI when all angular momentum are nonzero
   *
   *
   *  \param [in] pair1  bra shell pair data for shell1,shell2
   *  \param [in] pair2  ket shell pair data for shell3,shell4
   *  \param [in] shell1
   *  \param [in] shell2
   *  \param [in] shell3
   *  \param [in] shell4
   *  \param [in] FmT_2e Boys function between two shell pairs
   *  \param [in] LA
   *  \param [in] lA
   *  \param [in] LB
   *  \param [in] lB
   *  \param [in] LC
   *  \param [in] lC
   *  \param [in] LD
   *  \param [in] lD
   *
   *  \return ERI of two shell pairs ( shell1 shell2 | shell3 shell4 )
   */
  //----------------------------------------------------//
  // two-e horizontal recursion from (ab|cd) to (a0|cd) //
  // (ab|cd)=(a+1,b-1|cd)+(A-B)*(a,b-1|cd)              //
  //----------------------------------------------------//
   
  double RealGTOIntEngine::twoehRRabcd( 
     libint2::ShellPair &pair1 ,libint2::ShellPair &pair2 ,
     libint2::Shell &shell1, libint2::Shell &shell2,
     libint2::Shell &shell3, libint2::Shell &shell4,
     std::vector<std::vector<std::vector<double>>> &FmT_2e,
     int LA,int *lA,int LB,int *lB,int LC,int *lC,int LD,int *lD) {

     double tmpVal = 0.0, tmpVal1=0.0;

  // iWork is used to indicate which Cartesian angular momentum we are reducing (x,y,z)

    int iWork;
    int totalL = LA + LB + LC + LD;


    if(totalL==0) {   // (SS||SS) type 

      return twoeSSSS0(pair1,pair2,shell1,shell2,shell3,shell4);

      } else { // when totalL! = 0

        if( LB>=1 ) {

          int lAp1[3],lBm1[3];
          for( iWork=0 ; iWork<3 ; iWork++ ){
            lAp1[iWork]=lA[iWork];     
            lBm1[iWork]=lB[iWork];     
          } // for iWork

          if (lB[0]>0)      iWork=0;   
          else if (lB[1]>0) iWork=1;
          else if (lB[2]>0) iWork=2;
          lAp1[iWork]+=1;
          lBm1[iWork]-=1;

          tmpVal += twoehRRabcd(pair1,pair2,shell1,shell2,shell3,shell4,
                      FmT_2e, LA+1,lAp1,LB-1,lBm1,LC,lC,LD,lD);

          if ( std::abs(pair1.AB[iWork]) > 1.0e-15 ) {
            tmpVal += pair1.AB[iWork]*twoehRRabcd(pair1,pair2,shell1,shell2,shell3,
                        shell4,FmT_2e, LA,lA,LB-1,lBm1,LC,lC,LD,lD);
          } 

        } else if ( LB == 0 ) {
          tmpVal = twoehRRa0cd(pair1,pair2,shell1,shell2,shell3,shell4,FmT_2e,
                                 LA,lA,LC,lC,LD,lD);

        } // LB == 0
      }  // else ( LTOTAL != 0 )

      return tmpVal;

  }  // twoehRRabcd

  /**
   *  \brief horizontal recursion of ERI when all LB=0
   *
   *
   *  \param [in] pair1  bra shell pair data for shell1,shell2
   *  \param [in] pair2  ket shell pair data for shell3,shell4
   *  \param [in] shell1
   *  \param [in] shell2
   *  \param [in] shell3
   *  \param [in] shell4
   *  \param [in] FmT_2e Boys function between two shell pairs
   *  \param [in] LA
   *  \param [in] lA
   *  \param [in] LC
   *  \param [in] lC
   *  \param [in] LD
   *  \param [in] lD
   *
   *  \return ERI of two shell pairs ( shell1 shell2 | shell3 shell4 )
   */
  //----------------------------------------------------//
  // two-e horizontal recursion from (a0|cd) to (a0|c0) //
  // (a0|cd)=(a,0|c+1,d-1)+(C-D)*(a,0|c,d-1)            //
  //----------------------------------------------------//
   
  double RealGTOIntEngine::twoehRRa0cd(
    libint2::ShellPair &pair1, libint2::ShellPair &pair2,
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4, 
    std::vector<std::vector<std::vector<double>>> &FmT_2e,
    int LA,int *lA,int LC,int *lC,int LD,int *lD)  {

    double tmpVal=0.0;
    if(LD==0) {
      int pair1index=0, pair2index=0;
      // go into the vertical recursion
      for ( auto &pripair1 : pair1.primpairs ) {
        pair2index = 0; 
        for ( auto &pripair2 : pair2.primpairs ) {

          auto norm = 
                 shell1.contr[0].coeff[pripair1.p1]* 
                 shell2.contr[0].coeff[pripair1.p2]* 
                 shell3.contr[0].coeff[pripair2.p1]* 
                 shell4.contr[0].coeff[pripair2.p2];  

          tmpVal +=  norm * twoevRRa0c0( pripair1, pripair2,  
             FmT_2e[pair1index][pair2index], shell1,shell3, 0, LA,lA,LC,lC);

          pair2index++;
        } // for pripair2
        pair1index++;
      } // for pripair1
    } else { // if LD>0, go into horizontal recursion 

      int iWork;
      int lCp1[3],lDm1[3];  
     
      for( int iWork=0 ; iWork<3 ; iWork++ ){
        lCp1[iWork]=lC[iWork];
        lDm1[iWork]=lD[iWork];
      }
     
      if (lD[0]>0) iWork=0;
      else if (lD[1]>0) iWork=1;
      else if (lD[2]>0) iWork=2;
     
      lCp1[iWork]+=1;


    // when LD > 0

      lDm1[iWork] -=1 ;
      tmpVal = twoehRRa0cd(pair1, pair2, shell1, shell2, shell3, shell4, 
                             FmT_2e,LA,lA, LC+1,lCp1, LD-1,lDm1 ); 
      if ( std::abs(pair2.AB[iWork]) > 1.0e-15 ){
        tmpVal += pair2.AB[iWork] * twoehRRa0cd( pair1, pair2, shell1, shell2, 
               shell3, shell4, FmT_2e, LA,lA, LC,lC, LD-1,lDm1 );
      }

    } // else ( that means LD > 0 )
    return tmpVal;
  }  // twoehRRa0cd


  /**
   *  \brief vertical recursion of ERI when all LA, LC > 0
   *
   *
   *  \param [in] pripair1  primitive bra shell pair data for shell1,shell2
   *  \param [in] pripair2  primitive ket shell pair data for shell3,shell4
   *  \param [in] pair1index index of primitive shell pair among the contracted pair
   *  \param [in] pair2index index of primitive shell pair among the contracted pair
   *  \param [in] FmT_2e Boys function between two primitive shell pairs
   *  \param [in] shell1    
   *  \param [in] shell3    
   *  \param [in] m         order of auxiliary function
   *  \param [in] LA
   *  \param [in] lA
   *  \param [in] LC
   *  \param [in] lC
   *
   *  \return ERI of two primitive shell pairs ( shell1 shell2 | shell3 shell4 )
   */

//---------------------------------------------------------------//
// two-e vertical recursion from [a0|c0] to [a0|00]              //
// [a0|c0]^m = (Q-C)*[a0|c-1,0]^m                                //
//           + (W-Q)*[a0|c-1,0]^(m+1)                            //
//           + N(a)/(2*(zeta+eta))*[a-1,0|c-1,0]^(m+1)           //
//           + (N(c)-1)/(2*eta)*[a0|c-2,0]^m                     //
//           - (N(c)-1)/(2*eta)*zeta/(zeta+eta)*[a0|c-2,0]^(m+1) //
//---------------------------------------------------------------//
   
  double RealGTOIntEngine::twoevRRa0c0(
    libint2::ShellPair::PrimPairData &pripair1,
    libint2::ShellPair::PrimPairData &pripair2, 
    std::vector<double> &FmT_2epri, 
    libint2::Shell &shell1, libint2::Shell &shell3,
    int m, int LA, int *lA, int LC, int *lC ) {
 
 
    if(LC==0) return twoevRRa000( pripair1, pripair2, FmT_2epri,
                                  shell1, m, LA, lA );
 
    int lAm1[3],lCm1[3];  
 
    for ( int iWork=0 ; iWork<3 ; iWork++ ){
      lAm1[iWork]=lA[iWork];     
      lCm1[iWork]=lC[iWork];
    }
 
    double tmpVal=0.0;
    double W_iWork,Zeta;
    int iWork;
 
    if (lC[0]>0) iWork=0;
    else if (lC[1]>0) iWork=1;
    else if (lC[2]>0) iWork=2;

    lCm1[iWork]-=1;

    if ( std::abs(pripair2.P[iWork] - shell3.O[iWork]) > 1.0e-15 ) {
      tmpVal += ( pripair2.P[iWork] - shell3.O[iWork] ) * 
                twoevRRa0c0( pripair1, pripair2, FmT_2epri, 
                  shell1, shell3, m, LA,lA, LC-1,lCm1 );
    }  // if (Q-C) > 1.0e-15 


    Zeta = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma; 
/*
      for ( int mu = 0 ; mu<3 ; mu++ )
      W[mu] = (pripair1.P[mu]/pripair1.one_over_gamma 
              + pripair2.P[mu]/pripair2.one_over_gamma )/Zeta ; 
*/

    W_iWork = ( pripair1.P[iWork]/pripair1.one_over_gamma + 
               pripair2.P[iWork]/pripair2.one_over_gamma )/Zeta;

    if( std::abs( W_iWork- pripair2.P[iWork] ) > 1.0e-15 ) {
      tmpVal += ( W_iWork- pripair2.P[iWork] ) * twoevRRa0c0( pripair1, pripair2, 
                 FmT_2epri, shell1,shell3, m+1, LA,lA, LC-1,lCm1 );
    } // if( abs( W_iWork- Q[iWork] ) > 1.0e-15 )

    if (lA[iWork]>0) {

      lAm1[iWork] -= 1;
      tmpVal += (lAm1[iWork]+1) / (2.0*Zeta) * twoevRRa0c0( pripair1, pripair2, 
                 FmT_2epri, shell1, shell3, m+1, LA-1,lAm1, LC-1,lCm1 );
    } // if (lA[iWork]>0) 

    if ( lC[iWork]>=2 ){

      lCm1[iWork] -=1; // right now lCm1(iWork) = lC[iWork]-2 
      tmpVal += 0.5 * (lCm1[iWork]+1) * pripair2.one_over_gamma * 
        ( twoevRRa0c0( pripair1, pripair2, FmT_2epri, 
                       shell1, shell3, m, LA,lA, LC-2,lCm1 )

               - ( 1.0/pripair1.one_over_gamma )/Zeta 
          * twoevRRa0c0( pripair1, pripair2, FmT_2epri, 
                         shell1, shell3, m+1, LA,lA, LC-2,lCm1 ) );
    } // if ( lC[iWork]>=2 )

    return tmpVal;
  }  // twoevRRa0c0



  /**
   *  \brief vertical recursion of ERI when all LA > 0, all the others are 0
   *
   *
   *  \param [in] pripair1  primitive bra shell pair data for shell1,shell2
   *  \param [in] pripair2  primitive ket shell pair data for shell3,shell4
   *  \param [in] pair1index index of primitive shell pair among the contracted pair
   *  \param [in] pair2index index of primitive shell pair among the contracted pair
   *  \param [in] FmT_2e Boys function between two primitive shell pairs
   *  \param [in] shell1    
   *  \param [in] m         order of auxiliary function
   *  \param [in] LA
   *  \param [in] lA
   *
   *  \return ERI of two primitive shell pairs ( shell1 shell2 | shell3 shell4 )
   */
//---------------------------------------------------------------//
// two-e vertical recursion from [a0|00] to [00|00]              //
// [a0|00]^m = (P-A)*[a-1,0|00]^m                                //
//           + (W-P)*[a-1,0|00]^(m+1)                            //
//           + (N(a)-1)/(2*zeta)*[a-2,0|00]^m                    //
//           - (N(a)-1)/(2*zeta)*eta/(zeta+eta)*[a-2,0|00]^(m+1) //
//---------------------------------------------------------------//
   
  double RealGTOIntEngine::twoevRRa000(
    libint2::ShellPair::PrimPairData &pripair1,
    libint2::ShellPair::PrimPairData &pripair2, std::vector<double> &FmT_2epri,
    libint2::Shell &shell1, int m,int LA,int *lA ) {


    if(LA==0) {
      // calculate the (SS||SS) integral  	
      double expoT,Kab,Kcd,SSSS=0.0 ;
 
      // zeta+eta is expoT
      expoT = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma;  
 
      Kab = sqrt(2.0)* pow(M_PI,1.25) * pripair1.K; 
      Kcd = sqrt(2.0)* pow(M_PI,1.25) * pripair2.K;
 
      SSSS += Kab * Kcd * FmT_2epri[m] / sqrt(expoT); 
 
      return SSSS;
 
    } // if LA==0
 
    // here LA != 0
    double tmpVal=0.0,W[3],Zeta;
    int iWork;
    int lAm1[3];
 
    for( iWork=0 ; iWork<3 ; iWork++ ) lAm1[iWork]=lA[iWork];
 
    if (lA[0]>0) iWork=0;
    else if (lA[1]>0) iWork=1;
    else if (lA[2]>0) iWork=2;
 
    if( LA>=1 ) {
 
      lAm1[iWork]-=1;

      Zeta = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma;
      for ( int mu = 0 ; mu<3 ; mu++ )
        W[mu] = (pripair1.P[mu]/pripair1.one_over_gamma 
              + pripair2.P[mu]/pripair2.one_over_gamma )/Zeta ; 
 
      if ( std::abs( W[iWork]- pripair1.P[iWork] ) > 1.0e-15 ) {
 
        tmpVal += ( W[iWork]- pripair1.P[iWork] ) * twoevRRa000( pripair1, pripair2,
          FmT_2epri, shell1, m+1, LA-1,lAm1 );
 
      }  // if ( abs( W_iWork-P[iWork] )>1.0e-15 )

      if ( std::abs( pripair1.P[iWork] - shell1.O[iWork] )>1.0e-15 ) {  
        tmpVal+= ( pripair1.P[iWork] - shell1.O[iWork] ) * twoevRRa000( pripair1, 
          pripair2, FmT_2epri, shell1, m, LA-1,lAm1 );
      } // if ( abs( P[iWork]-A[iWork] )>1.0e-15 ) 

      if ( lA[iWork]>=2 ) {

        lAm1[iWork] -=1; // now lAm1[iWork] == lA[iWork]-2
        tmpVal += 0.5 * ( lAm1[iWork]+1 ) * pripair1.one_over_gamma  
                  *(twoevRRa000( pripair1, pripair2, FmT_2epri,
                    shell1, m, LA-2,lAm1 )
                - 1.0/(pripair2.one_over_gamma*Zeta) * twoevRRa000( pripair1, pripair2,
                  FmT_2epri, shell1, m+1, LA-2,lAm1 ) );
      } // if lA[iWork]>=2

    } // if( LA>=1 ) 

    return tmpVal;

  };  // twoevRRa000



  /**
   *  \brief vertical recursion of ERI when all the angular momentum are 0
   *
   *
   *  \param [in] pair1   bra shell pair data for shell1,shell2
   *  \param [in] pair2   ket shell pair data for shell3,shell4
   *  \param [in] shell1    
   *  \param [in] shell2    
   *  \param [in] shell3    
   *  \param [in] shell4    
   *
   *  \return ERI of two primitive shell pairs ( shell1 shell2 | shell3 shell4 )
   */
   
  double RealGTOIntEngine::twoeSSSS0(
    libint2::ShellPair &pair1, libint2::ShellPair &pair2,
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4 ) {

    // in this auxiliary integral, m=0 
    double norm,sqrPQ,PQ,expoT,T,FmT[1],Kab,Kcd,SSSS0=0.0 ;
    for ( auto &pripair1 : pair1.primpairs )
    for ( auto &pripair2 : pair2.primpairs ) {
       
      sqrPQ = 0.0;
      for(int m=0 ; m<3 ; m++ ) {
        PQ = ( pripair1.P[m]-pripair2.P[m] );
        sqrPQ += PQ*PQ;
      } 

      // zeta+eta is expoT
      expoT = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma;  

      // T= \rho * PQ^2
      T = sqrPQ/( pripair1.one_over_gamma + pripair2.one_over_gamma );

      // calculate F0(T)
      computeFmTTaylor( FmT, T, 0, 0 );

      Kab = sqrt(2.0)* pow(M_PI,1.25) * pripair1.K; 
      Kcd = sqrt(2.0)* pow(M_PI,1.25) * pripair2.K;

      norm = shell1.contr[0].coeff[pripair1.p1]* 
             shell2.contr[0].coeff[pripair1.p2]* 
             shell3.contr[0].coeff[pripair2.p1]* 
             shell4.contr[0].coeff[pripair2.p2];  

      SSSS0 += norm * Kab * Kcd * FmT[0] / sqrt(expoT); 

    } // for pripair2

    return SSSS0;
    
  } // twoeSSSS0

}  // namespace ChronusQ 

