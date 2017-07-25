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
#ifndef __INCLUDED_CQLINALG_BLAS1_HPP__
#define __INCLUDED_CQLINALG_BLAS1_HPP__

#include <cqlinalg/cqlinalg_config.hpp>

namespace ChronusQ {

  /**
   *  \brief Takes the inner product of two, possibly strided
   *  vectors
   *
   *  \f[
   *    a = \sum_i x_i^* y_i
   *  \f]
   *
   *  General to real, complex and mixed expressions.
   *
   *  \warning Real return involving complex vectors will
   *  return only the real part and discard the imaginary part.
   *
   *  Wraps BLAS functions. See
   *    http://www.netlib.org/lapack/lapack-3.1.1/html/ddot.f.html or
   *    http://www.netlib.org/lapack/lapack-3.1.1/html/zdotc.f.html for
   *  parameter documentation.
   */ 
  template <typename _F1, typename _F2, typename _F3>
  _F1 InnerProd(int N, _F2 *X, int INCX, _F3 *Y, int INCY);

  /**
   *  \brief Returns the euclidian norm of a vector
   *
   *  Wraps BLAS functions. See
   *    http://www.netlib.org/lapack/lapack-3.1.1/html/dnrm2.f.html or
   *    http://www.netlib.org/lapack/lapack-3.1.1/html/dznrm2.f.html for
   *  parameter documentation.
   *    
   */
  template <typename _F1, typename _F2>
  _F1 TwoNorm(int N, _F2 *X, int INCX);
  

}; // namespace ChronusQ


#endif

