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
#ifndef __INCLUDED_CQLINALG_SOLVE_HPP__
#define __INCLUDED_CQLINALG_SOLVE_HPP__

#include <cqlinalg/cqlinalg_config.hpp>

namespace ChronusQ {

  /**
   *  \brief Solves a linear system AX = B. Smart wrapper around
   *  DGESV or ZGESV depending on context.
   *
   *  Assumes preallocated IPIV
   *
   *  See http://www.netlib.org/lapack/lapack-3.1.1/html/dgesv.f.html or
   *      http://www.netlib.org/lapack/lapack-3.1.1/html/zgesv.f.html for
   *  parameter documentation.
   */ 
  template <typename _F>
  int LinSolve(int N, int NRHS, _F *A, int LDA, _F *B, 
    int LDB, int* IPIV);



  /**
   *  \brief Solves a linear system AX = B. Smart wrapper around
   *  DGESV or ZGESV depending on context.
   *
   *  Allocates memory internally through CQMemManager.
   *
   *  See http://www.netlib.org/lapack/lapack-3.1.1/html/dgesv.f.html or
   *      http://www.netlib.org/lapack/lapack-3.1.1/html/zgesv.f.html for
   *  parameter documentation.
   */ 
  template <typename _F>
  int LinSolve(int N, int NRHS, _F *A, int LDA, _F *B, 
    int LDB, CQMemManager &mem) {

    int* iPIV = mem.malloc<int>(N);

    LinSolve(N,NRHS,A,LDA,B,LDB,iPIV);

    mem.free(iPIV);

  };


  /**
   *  \brief Solves a linear system AX = B. Smart wrapper around
   *  DGESV or ZGESV depending on context.
   *
   *  Allocates on the stack
   *
   *  See http://www.netlib.org/lapack/lapack-3.1.1/html/dgesv.f.html or
   *      http://www.netlib.org/lapack/lapack-3.1.1/html/zgesv.f.html for
   *  parameter documentation.
   */ 
  template <typename _F>
  int LinSolve(int N, int NRHS, _F *A, int LDA, _F *B, 
    int LDB) {

    std::vector<int> iPIV(N,0);
    LinSolve(N,NRHS,A,LDA,B,LDB,&iPIV[0]);

  };

}; // namespace ChronusQ

#endif
