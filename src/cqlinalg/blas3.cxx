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
#include <cqlinalg/blas3.hpp>

namespace ChronusQ {

  template<>
  void Gemm(char TRANSA, char TRANSB, int M, int N, int K, double ALPHA,
    double *A, int LDA, double *B, int LDB, double BETA, double *C, int LDC){
#ifdef _CQ_MKL
    dgemm
#else
    dgemm_
#endif
    (&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);

  }; // GEMM (real,real,real)


  template<>
  void Gemm(char TRANSA, char TRANSB, int M, int N, int K, dcomplex ALPHA,
    dcomplex *A, int LDA, dcomplex *B, int LDB, dcomplex BETA, dcomplex *C, 
    int LDC){
#ifdef _CQ_MKL
    zgemm(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);
#else
    zgemm_(&TRANSA,&TRANSB,&M,&N,&K,reinterpret_cast<double*>(&ALPHA),
      reinterpret_cast<double*>(A),&LDA,reinterpret_cast<double*>(B),&LDB,
      reinterpret_cast<double*>(&BETA),reinterpret_cast<double*>(C),&LDC);
#endif

  }; // GEMM (complex,complex,complex)


  template<>
  void Gemm(char TRANSA, char TRANSB, int M, int N, int K, dcomplex ALPHA,
    double *A, int LDA, dcomplex *B, int LDB, dcomplex BETA, dcomplex *C, 
    int LDC){
#ifdef _CQ_MKL
    dzgemm(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,A,&LDA,B,&LDB,&BETA,C,&LDC);
#else
    assert( TRANSA == 'N' and (TRANSB == 'N' or TRANSB == 'C') );

    int COLS_B = (TRANSB == 'N') ? N : K;

    Eigen::Map<
      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
    > AMap(A,LDA,K);

    Eigen::Map<
      Eigen::Matrix<dcomplex,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
    > BMap(B,LDB,COLS_B), CMap(C,LDC,N);

    if(TRANSB == 'N')
      CMap.block(0,0,M,N).noalias() = 
        AMap.block(0,0,M,K).cast<dcomplex>() *
        BMap.block(0,0,K,N);
    else
      CMap.block(0,0,M,N).noalias() = 
        AMap.block(0,0,M,K).cast<dcomplex>() *
        BMap.block(0,0,N,K).adjoint();
#endif

  }; // GEMM (real,complex,complex)


}; // namespace ChronusQ
