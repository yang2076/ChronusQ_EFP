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


}; // namespace ChronusQ
