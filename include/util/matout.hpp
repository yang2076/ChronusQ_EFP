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
#ifndef __INCLUDED_UTIL_MATOUT_HPP__
#define __INCLUDED_UTIL_MATOUT_HPP__

#include <chronusq_sys.hpp>
#include <cxxapi/output.hpp>
#include <Eigen/Core>

namespace ChronusQ {

constexpr long double PRINT_SMALL = 1e-10;

template <typename T>
void prettyPrintSmartBase(std::ostream& out, const T* A, const size_t M, 
  const size_t N, const size_t LDA, const size_t colStride = 1, 
  const size_t list = 5, const size_t printWidth = 16) {

  int end,endLT;
  out << bannerTop;

  out << std::scientific << std::left << std::setprecision(8);
  for(size_t i = 0; i < N; i += list) {
    out << std::endl;
    end = list;
    out << std::setw(5) << " ";
    if((i + list) >= N) end = N - i;
    for(size_t k = i; k < i+end; k++) out << std::setw(printWidth) << k+1;
    out << std::endl;
    for(size_t j = 0; j < M; j++) {
      out << std::setw(5) << std::left << j+1;
      out << std::right;
      for(size_t n = i; n < i+end; n++) {
        T VAL = A[j*colStride + n*LDA];
        if(std::abs(VAL) > PRINT_SMALL) out << std::setw(printWidth) << VAL; 
        else if(std::isnan(VAL))        out << std::setw(printWidth) << "NAN";
        else if(std::isinf(VAL))        out << std::setw(printWidth) << "INF";
        else                            out << std::setw(printWidth) << 0.;
      }
      out << std::endl;
    };
  };
  out << bannerEnd << std::endl;

}; // prettyPrintSmartBase

template <typename T, 
          typename std::enable_if<
                     std::is_same<T,double>::value,int>::type = 0>
void prettyPrintSmart(std::ostream& out, std::string str, const T* A, 
  const size_t M, const size_t N, const size_t LDA, const size_t colStride = 1, 
  const size_t list = 5, const size_t printWidth = 16) {

  out.precision(10);
  out.fill(' ');
  out << std::endl << str + ": " << std::endl;

  prettyPrintSmartBase(out,A,M,N,LDA,colStride,list,printWidth);

}; // prettyPrintSmart (T = double)

template <typename T, 
          typename std::enable_if<
                     std::is_same<T,dcomplex>::value,int>::type = 0>
void prettyPrintSmart(std::ostream& out, std::string str, const T* A, 
  const size_t M, const size_t N, const size_t LDA, const size_t colStride = 1, 
  const size_t list = 5, const size_t printWidth = 16) {

  out.precision(10);
  out.fill(' ');

  out << std::endl << "Re[" << str + "]: " << std::endl;
  prettyPrintSmartBase(out,reinterpret_cast<double*>(A),M,N,LDA,2*colStride,
    list,printWidth);

  out << std::endl << "Im[" << str + "]: " << std::endl;
  prettyPrintSmartBase(out,reinterpret_cast<double*>(A)+1,M,N,LDA,2*colStride,
    list,printWidth);

}; // prettyPrintSmart (T = dcomplex)

}; // namespace ChronusQ

#endif
