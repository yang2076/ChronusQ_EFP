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
#include <cqlinalg/blas1.hpp>

namespace ChronusQ {

  template<>
  double InnerProd(int N, double *X, int INCX, double *Y, int INCY) {
    return 
#ifdef _CQ_MKL
      ddot
#else
      ddot_
#endif 
        (&N,X,&INCX,Y,&INCY);
  }; // InnerProd real = (real,real)




  template<>
  dcomplex InnerProd(int N, double *X, int INCX, double *Y, int INCY) {
    return dcomplex(InnerProd<double>(N,X,INCX,Y,INCY));
  }; // InnerProd complex = (real,real)




  /**
   *  \warning Discards imaginary part of inner product
   */ 
  template<>
  double InnerProd(int N, dcomplex *X, int INCX, double *Y, int INCY) {
    INCX *= 2;
    return
#ifdef _CQ_MKL
      ddot
#else
      ddot_
#endif 
        (&N,reinterpret_cast<double*>(X),&INCX,Y,&INCY);
  }; // InnerProd real = (complex,real)




  template<>
  dcomplex InnerProd(int N, dcomplex *X, int INCX, double *Y, int INCY) {
    INCX *= 2;
    double re =
#ifdef _CQ_MKL
      ddot
#else
      ddot_
#endif 
        (&N,reinterpret_cast<double*>(X),&INCX,Y,&INCY);

    double im = -
#ifdef _CQ_MKL
      ddot
#else
      ddot_
#endif 
        (&N,reinterpret_cast<double*>(X) + 1,&INCX,Y,&INCY);

    return dcomplex(re,im);
  }; // InnerProd complex = (complex,real)




  template<>
  dcomplex InnerProd(int N, dcomplex *X, int INCX, dcomplex *Y, int INCY) {
#ifdef _CQ_MKL
      dcomplex res;
      zdotc(&res,&N,X,&INCX,Y,&INCY);
      return res;
#else
      return zdotc_(&N,reinterpret_cast<double*>(X),&INCX,
        reinterpret_cast<double*>(Y),&INCY);
#endif
  }; // InnerProd complex = (complex,complex)

  /**
   *  \warning Discards imaginary part of inner product
   */ 
  template<>
  double InnerProd(int N, dcomplex *X, int INCX, dcomplex *Y, int INCY) {
    return
#ifdef _CQ_MKL
      std::real(InnerProd<dcomplex>(N,X,INCX,Y,INCY));
#else
      __real__ zdotc_(&N,reinterpret_cast<double*>(X),&INCX,
        reinterpret_cast<double*>(Y),&INCY);
#endif
  }; // InnerProd real = (complex,complex)



  template<>
  double TwoNorm(int N, double *X, int INCX) {
    return
#ifdef _CQ_MKL
      dnrm2
#else
      dnrm2_
#endif
      (&N,X,&INCX);
  }; // TwoNorm double = (double)

  template<>
  double TwoNorm(int N, dcomplex *X, int INCX) {
    return
#ifdef _CQ_MKL
      dznrm2(&N,X,&INCX);
#else
      dznrm2_(&N,reinterpret_cast<double*>(X),&INCX);
#endif
      
  }; // TwoNorm double = (dcomplex)

};
