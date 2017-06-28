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
#include <cqlinalg/factorization.hpp>

namespace ChronusQ {

  // Cholesky specializations

  // Real wraps DPOTRF
  template <>
  int Cholesky(char UPLO, int N, double *A, int LDA) {
    int INFO;
    dpotrf_(&UPLO,&N,A,&LDA,&INFO);
    return INFO;
  }; // Cholesky (real)

  // Complex wraps ZPOTRF
  template <>
  int Cholesky(char UPLO, int N, dcomplex *A, int LDA) {
    int INFO;
    zpotrf_(&UPLO,&N,A,&LDA,&INFO);
    return INFO;
  }; // Cholesky (complex)


  // Real wraps DPOTRI
  template <>
  int CholeskyInv(char UPLO, int N, double *A, int LDA) {
    int INFO;
    dpotri_(&UPLO,&N,A,&LDA,&INFO);
    return INFO;
  }; // CholeskyInv (real)

  // Complex wraps ZPOTRI
  template <>
  int CholeskyInv(char UPLO, int N, dcomplex *A, int LDA) {
    int INFO;
    zpotri_(&UPLO,&N,A,&LDA,&INFO);
    return INFO;
  }; // CholeskyInv (complex)

}; // namespace ChronusQ
